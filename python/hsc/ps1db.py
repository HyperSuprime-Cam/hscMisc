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

from .buildAndCatalog import BuildAndCatalog

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

class BuildPS1(BuildAndCatalog):
    def filter(self, data):
        mag = numpy.ndarray((5, len(data)))
        err = numpy.ndarray((5, len(data)))

        quality = data.field("o_qualityFlag")
        extended = (quality & QUALITY['EXT']) == 0

        for i, f in enumerate(FILTERS):
            mean = data.field("o_" + f + "MeanPSFMag")
            meanErr = data.field("o_" + f + "MeanPSFMagErr")
            badMean = (mean == -999) | ((quality & QUALITY['GOOD']) == 0) | (mean > LIMITS_MEAN[f])

            stack = data.field("o_" + f + "StackPSFMag")
            stackErr = data.field("o_" + f + "StackPSFMagErr")
            badStack = (stack == -999) | ((quality & QUALITY['GOOD_STACK']) == 0) | (stack > LIMITS_STACK[f])

            mag[i,:] = numpy.where(badMean, numpy.where(extended | badStack, numpy.nan, stack), mean)
            err[i,:] = numpy.where(badMean, numpy.where(extended | badStack, numpy.nan, stackErr), meanErr)

        numBad = numpy.isnan(mag).sum(axis=0)
        isBad = numBad > 3
        isGood = numpy.logical_not(isBad)

        return dict([("id", data.o_objID[isGood]),
                     ("ra", data.o_ra[isGood]),
                     ("dec", data.o_dec[isGood])] +
                    [(f, mag[i,isGood]) for i,f in enumerate(FILTERS)] +
                    [(f + "_err", err[i,isGood]) for i,f in enumerate(FILTERS)]
                    )
