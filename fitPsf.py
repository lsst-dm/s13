#!/usr/bin/env python

import sys
import os
import numpy
import scipy.optimize
from matplotlib import pyplot

import lsst.afw.geom.ellipses
import lsst.afw.image
import lsst.afw.detection
import lsst.afw.display.utils
import lsst.meas.multifit
import lsst.shapelet

def read():
    psf = lsst.afw.detection.Psf.readFits(os.path.join(os.environ["S13_DIR"], "psf.fits"))
    return psf.computeKernelImage(), psf.computeShape()

def asSamples(image):
    bbox = image.getBBox(lsst.afw.image.PARENT)
    x, y = numpy.meshgrid(numpy.arange(bbox.getBeginX(), bbox.getEndX()),
                          numpy.arange(bbox.getBeginY(), bbox.getEndY()))
    w = image.getArray().flatten()
    s = numpy.zeros((w.size, 2), dtype=float)
    s[:,0] = x.flatten()
    s[:,1] = y.flatten()
    return s, w

def fitMixture(s, w, radii, nIterations=20, order=4):
    components = lsst.meas.multifit.Mixture2.ComponentList()
    mu = numpy.zeros(2, dtype=float)
    for radius in radii:
        sigma = numpy.identity(2, dtype=float) * radius**2
        components.append(lsst.meas.multifit.Mixture2.Component(1.0, mu, sigma))
    mixture = lsst.meas.multifit.Mixture2(components)
    for i in range(nIterations):
        mixture.updateEM(s, w)
    print mixture
    msf = lsst.shapelet.MultiShapeletFunction()
    for component in mixture:
        mu = component.getMu()
        sigma = component.getSigma()
        ellipse = lsst.afw.geom.ellipses.Ellipse(
            lsst.afw.geom.ellipses.Quadrupole(sigma[0,0], sigma[1,1], sigma[0,1]),
            lsst.afw.geom.Point2D(mu[0], mu[1])
            )
        sf = lsst.shapelet.ShapeletFunction(order, lsst.shapelet.HERMITE, ellipse)
        sf.getCoefficients()[0] = component.weight / lsst.shapelet.ShapeletFunction.FLUX_FACTOR
        msf.getElements().push_back(sf)
    return msf

def fitShapelets(image, msf):
    bbox = image.getBBox(lsst.afw.image.PARENT)
    x, y = numpy.meshgrid(numpy.arange(bbox.getBeginX(), bbox.getEndX(), dtype=float),
                          numpy.arange(bbox.getBeginY(), bbox.getEndY(), dtype=float))
    matrix = numpy.zeros((sum(lsst.shapelet.computeSize(element.getOrder())
                              for element in msf.getElements()),
                          bbox.getArea()),
                         dtype=float).transpose()
    offset = 0
    for element in msf.getElements():
        builder = lsst.shapelet.ModelBuilderD(x.flatten(), y.flatten())
        size = lsst.shapelet.computeSize(element.getOrder())

        builder.addModelMatrix(element.getOrder(), matrix[:,offset:offset+size])
        offset += size
    data = image.getArray().flatten()
    coefficients, _, _, _ = numpy.linalg.lstsq(matrix, data)
    offset = 0
    for element in msf.getElements():
        size = lsst.shapelet.computeSize(element.getOrder())
        element.getCoefficients()[:] = coefficients[offset:offset+size]
        offset += size

def displayResiduals(image, msf):
    data = image.clone()
    factor = data.getArray().max()
    data /= factor
    model = data.clone()
    model.set(0.0)
    msf.evaluate().addToImage(model)
    model /= factor
    residuals = data.clone()
    residuals -= model
    print numpy.abs(residuals.getArray()).max()
    residuals *= 1000
    mosaic = lsst.afw.display.utils.Mosaic()
    mosaic.setMode("x")
    mosaic.append(data, "data")
    mosaic.append(model, "model")
    mosaic.append(residuals, "data-model")
    grid = mosaic.makeMosaic()
    lsst.afw.display.ds9.mtv(grid)

def fitFull(image, orders, radii=None, p0=None, penalty=0.01):
    bbox = image.getBBox(lsst.afw.image.PARENT)
    x, y = numpy.meshgrid(numpy.arange(bbox.getBeginX(), bbox.getEndX(), dtype=float),
                          numpy.arange(bbox.getBeginY(), bbox.getEndY(), dtype=float))
    basisOffsets = numpy.cumsum([0] + [lsst.shapelet.computeSize(order) for order in orders])
    matrix = numpy.zeros((basisOffsets[-1], bbox.getArea()),
                         dtype=float).transpose()
    builder = lsst.shapelet.ModelBuilderD(x.flatten(), y.flatten())
    data = image.getArray().flatten()

    if penalty is not None:
        penalty = numpy.identity(basisOffsets[-1], dtype=float) * penalty**2
        for i, order in enumerate(orders):
            penalty[basisOffsets[i], basisOffsets[i]] = 0.0

    if p0 is None:
        nRadii = len(radii)
        unitcircle = lsst.afw.geom.ellipses.Ellipse(
            lsst.afw.geom.ellipses.SeparableConformalShearLogTraceRadius(),
            lsst.afw.geom.Point2D()
            )
        p0 = numpy.zeros(nRadii*5, dtype=float)
        for i, radius in enumerate(radii):
            ellipse = lsst.afw.geom.ellipses.Ellipse(unitcircle)
            ellipse.scale(radius)
            p0[i*5:(i+1)*5] = ellipse.getParameterVector()

    def solve(p, *args):
        matrix[:] = 0.0
        for i, order in enumerate(orders):
            ellipse.setParameterVector(p[i*5:(i+1)*5])
            builder.update(ellipse)
            builder.addModelMatrix(order, matrix[:,basisOffsets[i]:basisOffsets[i+1]])
            f = numpy.dot(matrix.transpose(), matrix)
            g = numpy.dot(matrix.transpose(), data)
            if penalty is not None:
                f += penalty
        c, _, _, _  = numpy.linalg.lstsq(f, g)
        return c

    def func(p, *args):
        c = solve(p)
        r = numpy.dot(matrix, c)
        r -= data
        return r

    p1, flags = scipy.optimize.leastsq(func, p0, maxfev=10000)
    c1 = solve(p1)
    msf = lsst.shapelet.MultiShapeletFunction()
    for i, order in enumerate(orders):
        ellipse.setParameterVector(p1[i*5:(i+1)*5])
        sf = lsst.shapelet.ShapeletFunction(order, lsst.shapelet.HERMITE, ellipse,
                                            c1[basisOffsets[i]:basisOffsets[i+1]])
        msf.getElements().append(sf)
    return msf, func(p0)

def showMSF(msf):
    for element in msf.getElements():
        print element.getEllipse().getCore().getTraceRadius(), element.getCoefficients()
