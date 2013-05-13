#!/usr/bin/env python

import sys
import numpy
import pyfits

FILTERS = "ugriz"

def subset(inName, outName):
    inFile = pyfits.open(inName)
    inData = inFile[1].data

    schema = pyfits.ColDefs([pyfits.Column(name="id", format="K"),
                             pyfits.Column(name="thing_id", format="K"),
                             pyfits.Column(name="ra", format="D"),
                             pyfits.Column(name="dec", format="D"),
                             pyfits.Column(name="starnotgal", format="L"),
                             ] +
                            [pyfits.Column(name=name, format="E") for name in FILTERS] +
                            [pyfits.Column(name=name + "_err", format="E") for name in FILTERS]
                            )

    # From http://www.sdss3.org/dr9/algorithms/bitmask_calib_status.php
    # 0x0001 PHOTOMETRIC
    # 0x0002 OVERLAP
    # 0x0004 EXTRAP_CLEAR
    # 0x0008 EXTRAP_CLOUDY
    # 0x0010 DISJOINT
    # 0x0020 INCREMENT
    # 0x0040 RESERVED
    # 0x0080 RESERVED
    # 0x0100 PT_CLEAR
    # 0x0200 PT_CLOUDY
    # 0x0400 DEFAULT
    # 0x0800 NO_UBERCAL
    # Higher bits???
    indices = numpy.where(numpy.all(numpy.bitwise_and(inData.field("CALIB_STATUS"), 0x0001) > 0, axis=1))[0]
    if len(indices) == 0:
        print "No outputs for %s" % outName
        inFile.close()
        return

    outHdu = pyfits.new_table(schema, nrows=len(indices))
    outData = outHdu.data

    outData.id = inData.field("ID")[indices]
    outData.thing_id = inData.field("THING_ID")[indices]
    outData.ra = inData.field("RA")[indices]
    outData.dec = inData.field("DEC")[indices]
    for i, f in enumerate(FILTERS):
        # Some of the below (e.g., "mean") are functions in the pyfits.FITS_rec class,
        # so we need to access them differently than just grabbing an attribute.
        mean = outData.field(f)
        err = outData.field(f + "_err")

        mean[:] = inData.field(f.upper())[indices]
        err[:] = inData.field(f.upper() + "_ERR")[indices]

    outHdu.writeto(outName, clobber=True)
    print "Wrote %s" % outName
    inFile.close()

if __name__ == "__main__":
    for inName in sys.argv[1:]:
        subset(inName, "subset_" + inName)
