#!/usr/bin/env python

from argparse import ArgumentParser
import asciitable
import numpy
import pyfits
from collections import OrderedDict

FILTERS = "uBVgriz"
MAPPING = OrderedDict()
MAPPING['ID'] = pyfits.Column(name="id", format="K")
MAPPING['RA'] = pyfits.Column(name="ra", format="D")
MAPPING['DEC'] = pyfits.Column(name="dec", format="D")
MAPPING['star'] = pyfits.Column(name="starnotgal", format="L")
for f in FILTERS:
    MAPPING[f + '_mag'] = pyfits.Column(name=f, format="E")
    MAPPING[f + '_mag_err'] = pyfits.Column(name=f + "_err", format="E")


def convertCosmos(inName, outName):
    inFile = open(inName, "r")
    table = asciitable.read(inFile, Reader=asciitable.FixedWidthTwoLine, delimiter='|', header_start=0,
                            data_start=4, data_end=-1)

    schema = pyfits.ColDefs([column for column in MAPPING.values()])
    outHdu = pyfits.new_table(schema, nrows=len(table))
    outData = outHdu.data

    for name, column in MAPPING.items():
        outData.field(column.name)[:] = table.field(name)

    for f in FILTERS:
        mag = outData.field(f)
        err = outData.field(f + "_err")
        indices = numpy.where(numpy.logical_or(mag < 0, mag > 50))
        mag[indices] = numpy.NAN
        err[indices] = numpy.NAN

    outHdu.writeto(outName, clobber=True)
    print "Wrote %s" % outName
    print "To create an astrometry.net catalogue, execute:"
    outBase = outName.replace(".fits", "")
    print "build-index -i %s -o %s_and_0.fits -I 77770 -P0 -n 100 -S r -L 20 -E -M -j 0.4" % (inName, outBase)
    for i in range(1, 5):
        print "build-index -1 %s_and_0.fits -o %s_and_%d.fits -I 7777%d -P%d -n 100 -S r -L 10 -E -M -j 0.4 &" % (outBase, outBase, i, i, i)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("inFile", help="Name of input file")
    parser.add_argument("outFile", help="Name of output file")
    args = parser.parse_args()
    convertCosmos(args.inFile, args.outFile)
