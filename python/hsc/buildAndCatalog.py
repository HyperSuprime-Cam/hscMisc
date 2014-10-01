import re
import os
import glob
import multiprocessing
from argparse import ArgumentParser

import numpy
import pyfits


__all__ = ["BuildAndCatalog",]

# Support pickling of instance methods
import copy_reg, types
def unpickleInstanceMethod(obj, name):
    """Unpickle an instance method

    This has to be a named function rather than a lambda because
    pickle needs to find it.
    """
    return getattr(obj, name)
def pickleInstanceMethod(method):
    """Pickle an instance method

    The instance method is divided into the object and the
    method name.
    """
    obj = method.__self__
    name = method.__name__
    return unpickleInstanceMethod, (obj, name)
copy_reg.pickle(types.MethodType, pickleInstanceMethod)

class _Caller(object):
    def __init__(self, func):
        self.func = func
    def __call__(self, args):
        return self.func(*args[0], **args[1])
class MapFunc(object):
    """Class to wrap the map() function with optional threading"""
    def __init__(self, threads, func, maxtasksperchild=1):
        self.func = func
        self.args = []
        if threads > 1:
            self.pool = multiprocessing.Pool(threads, maxtasksperchild=maxtasksperchild)
            self.map = self.pool.map
        else:
            self.map = map
    def add(self, *args, **kwargs):
        self.args.append((args, kwargs))
    def run(self):
        caller = _Caller(self.func)
        return self.map(caller, self.args)

def system(command):
    print command
    os.system(command)


class BuildAndCatalog(object):
    def __init__(self, inputList, outputRoot, threads=0, nside=32):
        """Constructor

        The schema needs to be set appropriately for the output.  It is a
        dict with keys being the column names and values being the FITS
        column type.

        The build arguments needs to be set appropriately for the output.
        These are arguments to build-astrometry-index.
        """
        self.inputList = inputList
        self.outputRoot = outputRoot
        self.threads = threads
        self.nside = nside
        filters = "grizy"
        self.schema = dict([("id", "K"), ("ra", "D"), ("dec", "D")] + [(f, "E") for f in filters] +
                           [(f + "_err", "E") for f in filters])
        self.buildArgs = "-S i -L 20 -E -M -j 0.2 -n 100 -r 1"

    @classmethod
    def parse(cls):
        """Parse command-line arguments, returning a constructed BuildAndCatalog object."""
        parser = ArgumentParser()
        parser.add_argument("input", nargs="*", help="Input files")
        parser.add_argument("-j", dest="threads", type=int, default=0, help="Number of threads")
        parser.add_argument("-o", "--output", required=True, help="Output root name")
        parser.add_argument("-s", "--nside", type=int, default=32, help="HEALPix nside (power of 2)")
        args = parser.parse_args()
        return cls(args.input, args.output, threads=args.threads, nside=args.nside)

    def filter(self, data):
        """Filter the input data, returning the appropriate columns

        The subclass should define this to return a dict with keys
        being the column names and values being a numpy array of data.
        """
        raise NotImplementedError("Not implemented for base class")

    def convert(self, inName, outName):
        """Convert input data to the format to be processed by astrometry.net"""
        if os.path.exists(outName):
            print "Output file %s exists; not clobbering" % outName
            return
        inFile = pyfits.open(inName)
        inData = inFile[1].data
        print "Read %d rows from %s" % (len(inData), inName)

        # Filter the data and get the columns we want
        columns = self.filter(inData)

        if not "ra" in self.schema or not "dec" in self.schema:
            raise RuntimeError("Don't have 'ra' and 'dec' columns in schema")

        size = None
        for col in self.schema:
            if not col in columns:
                raise RuntimeError("Schema column %s was not present after filtering" % col)
            if size is None:
                size = len(columns[col])
            elif len(columns[col]) != size:
                raise RuntimeError("Size mismatch for column %s: %d vs %d" % (col, len(columns[col], size)))

        # Write it all out
        schema = pyfits.ColDefs([pyfits.Column(name=col, format=self.schema[col]) for col in self.schema])
        outHdu = pyfits.new_table(schema, nrows=size)
        outData = outHdu.data

        for col in self.schema:
            outData.field(col)[:] = columns[col]

        outHdu.writeto(outName, clobber=True)
        print "Wrote %d rows as %s" % (size, outName)
        inFile.close()

    def hpsplit(self, inputList):
        """Split the files into healpixes

        Only a single instance of this should be run on a catalog
        """
        out = "%s_hp_%%i.fits" % self.outputRoot
        args = "-r ra -d dec -n %d" % self.nside
        system("hpsplit -o " + out + " " + args + " " + " ".join(inputList))

    def generateIndexes(self, inName, index, healpix=None):
        """Generate astrometry.net indices

        Only a single instance of this should be run per input;
        inputs are usually divided into healpixes.
        """
        outName = "%s_and_%d" % (self.outputRoot, index)
        args = self.buildArgs[:] # Copy, so we're not overwriting when we append
        if healpix is not None:
            args += " -H %d" % healpix
        if self.nside is not None:
            args += " -s %d" % self.nside
        if not os.path.exists(outName + "_0.fits"):
            system("build-astrometry-index -i " + inName + " -o " + outName + "_0.fits -I " + str(index) + "0 -P 0 " + args)
        if not os.path.exists(outName + "_1.fits"):
            system("build-astrometry-index -1 " + outName + "_0.fits -o " + outName + "_1.fits -I " + str(index) + "1 -P 1 " + args)
        if not os.path.exists(outName + "_2.fits"):
            system("build-astrometry-index -1 " + outName + "_0.fits -o " + outName + "_2.fits -I " + str(index) + "2 -P 2 " + args)
        if False:
            # Don't need these: "-P 2  should work for images about 12 arcmin across" says build-astrometry-index
            system("build-astrometry-index -1 " + outName + "_0.fits -o " + outName + "_3.fits -I " + str(index) + "3 -P 3 " + args)
            system("build-astrometry-index -1 " + outName + "_0.fits -o " + outName + "_4.fits -I " + str(index) + "4 -P 4 " + args)

    def run(self):
        """Create astrometry.net indices

        Supports multiple 'threads' (though they're actually processes).
        """
        catList = []
        mapSubset = MapFunc(self.threads, self.convert)
        for i, inName in enumerate(self.inputList):
            catName = "%s_in_%d.fits" % (self.outputRoot, i)
            catList.append(catName)
            if not os.path.exists(catName):
                mapSubset.add(inName, catName)
        mapSubset.run()
        del mapSubset

        self.hpsplit(catList)

        mapIndexes = MapFunc(self.threads, self.generateIndexes)
        for inName in glob.glob("%s_hp_*.fits" % self.outputRoot):
            m = re.search(r"%s_hp_(\d+)\.fits" % self.outputRoot, inName)
            assert m, "Unable to match filename"
            healpix = int(m.group(1))
            mapIndexes.add(inName, healpix, healpix=healpix)
        mapIndexes.run()

    @classmethod
    def parseAndRun(cls):
        return cls.parse().run()
