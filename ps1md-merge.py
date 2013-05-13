#!/usr/bin/env python

import re
import sys
import numpy
import pyfits

# http://svn.pan-starrs.ifa.hawaii.edu/trac/ipp/wiki/MD.GR0#MD.V3Tessellation
SIZE = 6400 # Skycell size
BORDER = 320 # Skycell overlap
FILTERS = "grizy"
TEMPLATE = "i" # Template filter name

MASK = ["SOURCE.MASK." + s for s in ("BADPSF", "DEFECT", "CR_LIMIT", )]

def getFilter(inFits):
    return inFits[0].header["FPA.FILTER"].replace(".00000", "")

def getSkycell(filename):
    m = re.search(r"MD04\.V3\.skycell\.(...)\.sky", filename)
    assert m, "Unable to match"
    return m.groups()[0]

def emptyArray(num, dtype=numpy.float32):
    array = numpy.empty(num, dtype=numpy.float32)
    array[:] = numpy.NAN
    return array

def getStarGal(data):
    return (-2.5*numpy.log10(data.field("KRON_FLUX")) - data.field("PSF_INST_MAG")) > 0.0

class Data(object):
    def __init__(self, inFits):
        filterName = getFilter(inFits)
        data = inFits[1].data
        self.maskVal = sum(inFits[1].header[h] for h in MASK)
        self.size = (inFits[0].header["NAXIS1"], inFits[0].header["NAXIS2"])
        self.id = data.field("IPP_IDET").astype(numpy.int64)
        self.num = len(self.id)
        self.x = data.field("X_PSF").astype(numpy.float32)
        self.y = data.field("Y_PSF").astype(numpy.float32)
        self.ra = None
        self.dec = None
        self.flags = data.field("FLAGS").astype(numpy.int64)
        setattr(self, "mag_" + filterName, data.field("CAL_PSF_MAG").astype(numpy.float32))
        setattr(self, "mag_" + filterName + "_err", data.field("CAL_PSF_MAG_SIG").astype(numpy.float32))
        self.stargal = getStarGal(data)
        for f in set(FILTERS) - set(filterName):
            setattr(self, "mag_" + f, emptyArray(self.num))
            setattr(self, "mag_" + f + "_err", emptyArray(self.num))
        self.calculateRaDec(inFits)

    def calculateRaDec(self, inFits):
        header = inFits[0].header
        crval1 = header["CRVAL1"]
        crval2 = header["CRVAL2"]
        crpix1 = header["CRPIX1"]
        crpix2 = header["CRPIX2"]
        cdelt1 = header["CDELT1"]
        cdelt2 = header["CDELT2"]
        pc11 = header["PC001001"]
        pc12 = header["PC001002"]
        pc21 = header["PC002001"]
        pc22 = header["PC002002"]
        assert pc12 == 0.0 and pc21 == 0.0

        # Deprojection from psCoord.c
        xi = numpy.radians((self.x.astype(numpy.float64) - crpix1) * pc11 * cdelt1)
        eta = numpy.radians((self.y.astype(numpy.float64) - crpix2) * pc22 * cdelt2)
        r = numpy.hypot(xi, eta)
        rho = numpy.sqrt(1.0 + xi*xi + eta*eta)
        sinPhi = xi/r
        cosPhi = -eta/r
        sinTheta = 1.0/rho
        cosTheta = r/rho
        sinDp = numpy.sin(numpy.radians(crval2))
        cosDp = numpy.cos(numpy.radians(crval2))
        sinAlpha = cosTheta*sinPhi
        cosAlpha = cosTheta*cosPhi*sinDp + sinTheta*cosDp

        self.dec = numpy.degrees(numpy.arcsin(sinTheta*sinDp - cosTheta*cosPhi*cosDp))
        self.ra = numpy.degrees(numpy.arctan2(sinAlpha, cosAlpha)) + crval1

    def add(self, inFits):
        filterName = getFilter(inFits)
        data = inFits[1].data
        magArray = getattr(self, "mag_" + filterName)
        errArray = getattr(self, "mag_" + filterName + "_err")
        index = 0
        for i, mag, err, flag in zip(data.field("IPP_IDET").astype(numpy.int64),
                                           data.field("CAL_PSF_MAG").astype(numpy.float32),
                                           data.field("CAL_PSF_MAG_SIG").astype(numpy.float32),
                                           data.field("FLAGS").astype(numpy.int64),
                                     ):
            while index < self.num and self.id[index] < i:
                index += 1
            if index >= self.num:
                break
            if self.id[index] == i:
                magArray[index] = mag
                errArray[index] = err
                self.flags[index] = numpy.bitwise_or(self.flags[index], flag)

    def tossBad(self, border=BORDER):
        indices = numpy.where((self.x > border) & (self.x < self.size[0] - border) &
                              (self.y > border) & (self.y < self.size[1] - border) &
                              numpy.logical_not(numpy.bitwise_and(self.flags, self.maskVal)) &
                              numpy.logical_not(numpy.isnan(getattr(self, "mag_" + TEMPLATE))))
        for name in ["id", "x", "y", "ra", "dec", "flags",] + ["mag_" + f for f in FILTERS] + ["mag_" + f + "_err" for f in FILTERS]:
            setattr(self, name, getattr(self, name)[indices])
        self.num = len(self.id)


def merge(outName, inList):

    # Select sources detected in i-band
    skycellData = {}
    for inName in inList:
        inFile = pyfits.open(inName)
        filterName = getFilter(inFile)
        if filterName != TEMPLATE:
            continue
        skycell = getSkycell(inName)
        print inName, skycell, filterName
        skycellData[skycell] = Data(inFile)

    # Add other bands
    for inName in inList:
        inFile = pyfits.open(inName)
        filterName = getFilter(inFile)
        if filterName == TEMPLATE:
            continue
        skycell = getSkycell(inName)
        print inName, skycell, filterName
        skycellData[skycell].add(inFile)

    for data in skycellData.values():
        data.tossBad()

    # Merge and write
    num = sum(d.num for d in skycellData.values())
    schema = pyfits.ColDefs([pyfits.Column(name="id", format="K"),
                             pyfits.Column(name="skycell", format="J"),
                             pyfits.Column(name="ra", format="D"),
                             pyfits.Column(name="dec", format="D"),
                             pyfits.Column(name="starnotgal", format="L"),
                             ] +
                            [pyfits.Column(name=name, format="E") for name in FILTERS] +
                            [pyfits.Column(name=name + "_err", format="E") for name in FILTERS]
                            )

    outHdu = pyfits.new_table(schema, nrows=num)
    outData = outHdu.data
    outData.field("id")[:] = numpy.arange(num)

    maxId = 0
    for skycell, inData in skycellData.items():
        outIndices = numpy.where((outData.field("id") >= maxId) & (outData.field("id") < maxId + inData.num))
        maxId += inData.num

        fieldList = [("skycell", numpy.ones(inData.num, dtype=outData.field("skycell").dtype) * int(skycell)),
                     ("ra", inData.ra),
                     ("dec", inData.dec),
                     ("starnotgal", inData.stargal),
                     ] + \
                     [(f, getattr(inData, "mag_" + f)) for f in FILTERS] + \
                     [(f + "_err", getattr(inData, "mag_" + f + "_err")) for f in FILTERS]

        for field, value in fieldList:
            data = outData.field(field)
            data[outIndices] = value

    outHdu.writeto(outName, clobber=True)
    print "Wrote %s" % outName
    inFile.close()

if __name__ == "__main__":
    merge(sys.argv[1], sys.argv[2:])

    print "To build an astrometry.net catalogue, execute:"
    fileName = sys.argv[1]
    rootName =fileName.replace(".fits", "")
    print "build-index -i %s -o %s_and_0.fits -I 77770 -P0 -n 100 -S r -L 20 -E -M -j 0.4" % (fileName, rootName)
    for i in range(1, 5):
        print "build-index -1 %s_and_0.fits -o %s_and_%d.fits -I 77774 -P%d -n 100 -S r -L 20 -E -M -j 0.4" % (rootName, rootName, i, i)
