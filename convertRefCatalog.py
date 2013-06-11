#!/usr/bin/env python

import sys

import lsst.afw.table
import lsst.afw.geom
import lsst.afw.coord

class Keys(object):
    schema = lsst.afw.table.SimpleTable.makeMinimalSchema()
    sed = schema.addField("sed", type=str, doc="Filename of SED", size=32)
    mag_norm = schema.addField("mag.norm", type=float, doc="magnitude relative to SED (?)")
    redshift = schema.addField("redshift", type=float, doc="galaxy redshift")
    sindex = schema.addField("sindex", type=float, doc="Sersic index")
    ellipse_a = schema.addField("ellipse.a", type=float, doc="ellipse semi-major axis", units="arcsec")
    ellipse_b = schema.addField("ellipse.b", type=float, doc="ellipse semi-minor axis", units="arcsec")
    ellipse_theta = schema.addField("ellipse.theta", type=float, doc="ellipse position angle")
    mag_u = schema.addField("mag.u", type=float, doc="u-band magnitude")
    mag_g = schema.addField("mag.g", type=float, doc="g-band magnitude")
    mag_r = schema.addField("mag.r", type=float, doc="r-band magnitude")
    mag_i = schema.addField("mag.i", type=float, doc="i-band magnitude")
    mag_z = schema.addField("mag.z", type=float, doc="z-band magnitude")
    mag_y = schema.addField("mag.y", type=float, doc="y-band magnitude")

def main(inputPath, outputPath):
    catalog = lsst.afw.table.SimpleCatalog(Keys.schema)
    with open(inputPath, 'r') as inputFile:
        for n, line in enumerate(inputFile):
            if n == 0: continue
            data = line.split(",")
            record = catalog.addNew()
            record.setId(int(data[0]))
            coord = lsst.afw.coord.IcrsCoord(float(data[1])*lsst.afw.geom.degrees,
                                             float(data[2])*lsst.afw.geom.degrees)
            record.setCoord(coord)
            record.setString(Keys.sed, data[3])
            record.setD(Keys.mag_norm, float(data[4]))
            record.setD(Keys.redshift, float(data[5]))
            record.setD(Keys.sindex, float(data[6]))
            record.setD(Keys.ellipse_a, float(data[7]))
            record.setD(Keys.ellipse_b, float(data[8]))
            record.setD(Keys.ellipse_theta, float(data[9]))
            record.setD(Keys.mag_u, float(data[10]))
            record.setD(Keys.mag_g, float(data[11]))
            record.setD(Keys.mag_r, float(data[12]))
            record.setD(Keys.mag_i, float(data[13]))
            record.setD(Keys.mag_z, float(data[14]))
            record.setD(Keys.mag_y, float(data[15]))
    catalog.writeFits(outputPath)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
