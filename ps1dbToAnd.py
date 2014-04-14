#!/usr/bin/env python

"""
Convert results from Pan-STARRS "PSPS" database query to astrometry.net format

Here's the sort of "SELECT" statement the code is intended to handle:

SELECT o.objID AS o_objID, o.ra AS o_ra, o.dec AS o_dec, o.gMeanPSFMag AS o_gMeanPSFMag, o.gMeanPSFMagErr AS o_gMeanPSFMagErr, o.rMeanPSFMag AS o_rMeanPSFMag, o.rMeanPSFMagErr AS o_rMeanPSFMagErr, o.iMeanPSFMag AS o_iMeanPSFMag, o.iMeanPSFMagErr AS o_iMeanPSFMagErr, o.zMeanPSFMag AS o_zMeanPSFMag, o.zMeanPSFMagErr AS o_zMeanPSFMagErr, o.yMeanPSFMag AS o_yMeanPSFMag, o.yMeanPSFMagErr AS o_yMeanPSFMagErr, o.gStackPSFMag AS o_gStackPSFMag, o.gStackPSFMagErr AS o_gStackPSFMagErr, o.gFlags AS o_gFlags, o.rStackPSFMag AS o_rStackPSFMag, o.rStackPSFMagErr AS o_rStackPSFMagErr, o.rFlags AS o_rFlags, o.iStackPSFMag AS o_iStackPSFMag, o.iStackPSFMagErr AS o_iStackPSFMagErr, o.iFlags AS o_iFlags, o.zStackPSFMag AS o_zStackPSFMag, o.zStackPSFMagErr AS o_zStackPSFMagErr, o.zFlags AS o_zFlags, o.yStackPSFMag AS o_yStackPSFMag, o.yStackPSFMagErr AS o_yStackPSFMagErr, o.yFlags AS o_yFlags, o.objInfoFlag AS o_objInfoFlag, o.qualityFlag AS o_qualityFlag
INTO mydb.whatever
FROM object AS o

Add the following to restrict to the appropriate HSC survey area:

* hsc_fall:
WHERE ((o.ra >= 22*15 - 0.5 OR o.ra <= (2 + 40/60)*15 + 0.5) AND o.dec >= -1.5 AND o.dec <= 7.5) OR
(o.ra >= (1 + 50/60)*15 - 0.5 AND o.ra <= (2 + 40/60)*15 + 0.5 AND o.dec >= -7.5 AND o.dec <= -0.5)

* hsc_spring:
WHERE (o.ra >= 8.5*15 - 0.5 AND o.ra <= 15*15 + 0.5 AND o.dec >= -2.5 AND o.dec <= 5.5)

* hsc_north:
WHERE (o.ra >= (13 + 20/60)*15 - 0.5 AND o.ra <= (16 + 40/60)*15 + 0.5 AND o.dec >= 42 AND o.dec <= 44.5)

* elaisn1:
JOIN fgetNearbyObjEq((16+10/60)*15, 54, 2.5) AS cone ON cone.objid = o.objID


The hsc_fall and hsc_spring areas are too large to retrieve everything in one go, so need to divide into
parts.  This is simply accomplished by adding a WHERE clause on the modulus of the object id (
e.g., "o.objID % 3 == 1").  For a 10 GB limit on PV1.2, it was necessary to divide hsc_fall and hsc_spring
into three.
"""

import re
import os
import glob
import multiprocessing
from argparse import ArgumentParser

import numpy
import pyfits

FILTERS = "grizy"

QUALITY = {'EXT': 0x0001, # extended in our data (eg, PS)
           'EXT_ALT': 0x0002, # extended in external data (eg, 2MASS)
           'GOOD': 0x0004, # good-quality measurement in our data (eg,PS)
           'GOOD_ALT': 0x0008, # good-quality measurement in external data (eg, 2MASS)
           'GOOD_STACK': 0x0010, # good-quality object in the stack (> 1 good stack)
           'SUSPECT_STACK': 0x0020, # suspect object in the stack (> 1 good or suspect stack, less tham 2 good)
           'BAD_STACK': 0x0040, # good-quality object in the stack (> 1 good stack)
           }

LIMITS_MEAN = {'g': 20.5,
               'r': 21.0,
               'i': 21.0,
               'z': 20.0,
               'y': 19.0,
               }
LIMITS_STACK = {'g': 21.5,
                'r': 22.0,
                'i': 22.0,
                'z': 21.0,
                'y': 20.0,
                }

def convert(inName, outName):
    if os.path.exists(outName):
        print "Output file %s exists; not clobbering" % outName
        return
    inFile = pyfits.open(inName)
    inData = inFile[1].data
    print "Read %d rows from %s" % (len(inData), inName)

    # Throw out the bad stuff
    mag = numpy.ndarray((5, len(inData)))
    err = numpy.ndarray((5, len(inData)))

    quality = inData.field("o_qualityFlag")
    extended = (quality & QUALITY['EXT']) == 0

    for i, f in enumerate(FILTERS):
        mean = inData.field("o_" + f + "MeanPSFMag")
        meanErr = inData.field("o_" + f + "MeanPSFMagErr")
        badMean = (mean == -999) | ((quality & QUALITY['GOOD']) == 0) | (mean > LIMITS_MEAN[f])

        stack = inData.field("o_" + f + "StackPSFMag")
        stackErr = inData.field("o_" + f + "StackPSFMagErr")
        badStack = (stack == -999) | ((quality & QUALITY['GOOD_STACK']) == 0) | (stack > LIMITS_STACK[f])

        mag[i,:] = numpy.where(badMean, numpy.where(extended | badStack, numpy.nan, stack), mean)
        err[i,:] = numpy.where(badMean, numpy.where(extended | badStack, numpy.nan, stackErr), meanErr)

    numBad = numpy.isnan(mag).sum(axis=0)
    isBad = numBad > 3
    isGood = numpy.logical_not(isBad)

    ident = inData.o_objID[isGood]
    ra = inData.o_ra[isGood]
    dec = inData.o_dec[isGood]
    mag = mag[:,isGood]
    err = err[:,isGood]

    # Write it all out
    schema = pyfits.ColDefs([pyfits.Column(name="id", format="K"),
                             pyfits.Column(name="ra", format="D"),
                             pyfits.Column(name="dec", format="D")] +
                            [pyfits.Column(name=name, format="E") for name in FILTERS] +
                            [pyfits.Column(name=name + "_err", format="E") for name in FILTERS]
                            )

    outHdu = pyfits.new_table(schema, nrows=len(ident))
    outData = outHdu.data

    outData.id = ident
    outData.ra = ra
    outData.dec = dec
    for i, f in enumerate(FILTERS):
        outValue = outData.field(f)
        outErr = outData.field(f + "_err")

        outValue[:] = mag[i,:]
        outErr[:] = err[i,:]

    outHdu.writeto(outName, clobber=True)
    print "Wrote %d rows as %s" % (len(ident), outName)
    inFile.close()

def system(command):
    print command
    os.system(command)


def hpsplit(inList, outroot, nside=16):
    args = "-r ra -d dec -n %d" % nside
    system("hpsplit -o " + outroot + " " + args + " " + " ".join(inList))

def generateIndexes(inName, outName, index, healpix=None, nside=None):
    args = "-S r -L 20 -E -M -j 0.4 -n 100"
    if healpix is not None:
        args += " -H %d" % healpix
    if nside is not None:
        args += " -s %d" % nside
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

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("input", nargs="*", help="Input files")
    parser.add_argument("-j", dest="threads", type=int, default=0, help="Number of threads")
    parser.add_argument("-o", "--output", required=True, help="Output root name")
    parser.add_argument("-s", "--nside", default=32, help="HEALPix nside (power of 2)")
    args = parser.parse_args()

    ps1List = []
    mapSubset = MapFunc(args.threads, convert)
    for i, inName in enumerate(args.input):
        ps1Name = "%s_in_%d.fits" % (args.output, i)
        ps1List.append(ps1Name)
        if not os.path.exists(ps1Name):
            mapSubset.add(inName, ps1Name)
    mapSubset.run()
    del mapSubset

    hpsplit(ps1List, "%s_hp_%%i.fits" % args.output, nside=args.nside)

    mapIndexes = MapFunc(args.threads, generateIndexes)
    for inName in glob.glob("%s_hp_*.fits" % args.output):
        m = re.search(r"%s_hp_(\d+)\.fits" % args.output, inName)
        assert m, "Unable to match filename"
        healpix = int(m.group(1))
        andName = "%s_and_%d" % (args.output, healpix)
        mapIndexes.add(inName, andName, healpix, nside=args.nside, healpix=healpix)
    mapIndexes.run()
