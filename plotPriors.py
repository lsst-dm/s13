#!/usr/bin/env python

import sys
import os
import cPickle
import numpy
from matplotlib import pyplot
from matplotlib.colors import LogNorm
from matplotlib.colorbar import Colorbar

import lsst.afw.geom as geom
import lsst.afw.geom.ellipses as el
import lsst.meas.multifit

def convert(subset):
    dtype = numpy.dtype([("id", (str, 64)), ("sindex", float), ("a", float), ("b", float), ("theta", float),
                         ("ixx", float), ("iyy", float), ("ixy", float),
                         ("e1", float), ("e2", float), ("r", float)])
    inFile = os.path.join(os.environ["S13_DATA_DIR"], "reference_cats_dithered",
                          "%s_ellipticity.dat" % subset)
    m = numpy.loadtxt(inFile, skiprows=1, dtype=dtype, delimiter=',')

    ellipticity = numpy.zeros(m.shape, dtype=float)
    radius = numpy.zeros(m.shape, dtype=float)

    for i, r in enumerate(m):
        axes = el.Axes(r["a"], r["b"], (r["theta"] * geom.degrees).asRadians(), True)
        separable = el.SeparableConformalShearLogTraceRadius(axes)
        ellipticity[i] = separable.getEllipticity().getE()
        radius[i] = separable.getRadius()

    outFile = os.path.join(os.environ["S13_DATA_DIR"], "reference_cats_dithered", "%s.p" % subset)
    with open(outFile, 'w') as outStream:
        cPickle.dump((ellipticity, radius), outStream, protocol=2)

def read(subset):
    inFile = os.path.join(os.environ["S13_DATA_DIR"], "reference_cats_dithered", "%s.p" % subset)
    with open(inFile, 'r') as inStream:
        ellipticity, radius = cPickle.load(inStream)
    return ellipticity, radius

def transform3d(ellipticity, radius, m=4):
    n = ellipticity.size
    data = numpy.zeros((m*n, 3), dtype=float)
    theta = numpy.random.rand(n) * numpy.pi
    for i in range(m):
        data[n*i:n*(i+1),0] = ellipticity * numpy.cos(2*(theta + i*numpy.pi/m))
        data[n*i:n*(i+1),1] = ellipticity * numpy.sin(2*(theta + i*numpy.pi/m))
        data[n*i:n*(i+1),2] = radius
    return data

def plotDensity(data, ix, iy, xlabel, ylabel, mixture=None, bins=50, weights=None, cmap=pyplot.cm.Greys):
    h2d, xEdges, yEdges = numpy.histogram2d(data[:,ix], data[:,iy], normed=True, bins=bins, weights=weights)
    xCenters = 0.5*(xEdges[1:] + xEdges[:-1])
    yCenters = 0.5*(yEdges[1:] + yEdges[:-1])
    norm = LogNorm(clip=True)
    levels = [0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8]
    fig = pyplot.figure("%s vs %s" % (ylabel, xlabel))
    ax2d = fig.add_subplot(2,2,3)
    image = ax2d.imshow(h2d.transpose(), interpolation='nearest', origin='lower', cmap=cmap, alpha=0.8,
                        extent=(xEdges[0], xEdges[-1], yEdges[0], yEdges[-01]), norm=norm, aspect='auto')
    contours = ax2d.contour(xCenters, yCenters, h2d.transpose(), colors='k', norm=norm, levels=levels)
    ax2d.set_ylabel(ylabel)
    ax2d.set_xlabel(xlabel)
    ax1y = fig.add_subplot(2,2,4, sharey=ax2d)
    ax1y.hist(data[:,iy], bins=bins, weights=weights, orientation='horizontal', linewidth=0,
              facecolor='k', alpha=0.5, normed=True)
    ax1y.get_yaxis().set_visible(False)
    ax1x = fig.add_subplot(2,2,1, sharex=ax2d)
    ax1x.hist(data[:,ix], bins=bins, weights=weights, orientation='vertical', linewidth=0,
              facecolor='k', alpha=0.5, normed=True)
    ax1x.get_xaxis().set_visible(False)
    bbox1x = ax1x.get_position()
    bbox1y = ax1y.get_position()
    rect = (bbox1y.x0, bbox1x.y0, 0.2*(bbox1y.x1 - bbox1y.x0), bbox1x.y1 - bbox1x.y0)
    axcb = fig.add_axes(rect)
    cb = fig.colorbar(image, cax=axcb)
    cb.add_lines(contours)
    if mixture is not None:
        mix2d = mixture.project(ix, iy)
        xr = numpy.linspace(xEdges[0], xEdges[-1], 500)
        yr = numpy.linspace(yEdges[0], yEdges[-1], 500)
        xg, yg = numpy.meshgrid(xr, yr)
        xy = numpy.array([xg.flatten(), yg.flatten()]).transpose().copy()
        pxy = numpy.zeros(xy.shape[0], dtype=float)
        mix2d.evaluate(xy, pxy)
        pg = pxy.reshape(xg.shape)
        ax2d.contour(xg, yg, pg, colors='r', norm=norm, levels=levels)
        mix1x = mixture.project(ix)
        px = numpy.zeros(xr.shape, dtype=float)
        mix1x.evaluate(xr.reshape(-1,1), px)
        ax1x.plot(xr, px, 'r')
        mix1y = mixture.project(iy)
        py = numpy.zeros(yr.shape, dtype=float)
        mix1y.evaluate(yr.reshape(-1,1), py)
        ax1y.plot(py, yr, 'r')
    return ax2d, ax1x, ax1y

if __name__ == "__main__":
    if sys.argv[1] == "convert":
        convert(sys.argv[2])
    elif sys.argv[1] == "plot":
        plot(sys.argv[2])
    elif sys.argv[1] == "fit":
        fit(sys.argv[2])

