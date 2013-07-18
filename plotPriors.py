#!/usr/bin/env python

import sys
import os
import cPickle
import numpy
from matplotlib import pyplot

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

def plot(subset):
    inFile = os.path.join(os.environ["S13_DATA_DIR"], "reference_cats_dithered", "%s.p" % subset)
    with open(inFile, 'r') as inStream:
        ellipticity, radius = cPickle.load(inStream)
    pyplot.figure()
    pyplot.hist(radius, bins=50)
    pyplot.figure()
    pyplot.hist(ellipticity, bins=50)
    pyplot.figure()
    pyplot.hexbin(ellipticity, radius, bins=50)
    pyplot.show()

def read(subset):
    inFile = os.path.join(os.environ["S13_DATA_DIR"], "reference_cats_dithered", "%s.p" % subset)
    with open(inFile, 'r') as inStream:
        ellipticity, radius = cPickle.load(inStream)
    return ellipticity, radius

def fit(ellipticity, radius, nComponents=10):
    data = transform(ellipticity, radius)
    from sklearn import mixture
    clf = mixture.GMM(n_components=nComponents, cvtype='full')
    clf.fit(data)
    return data, clf

def plotFit(data, clf):
    cmap = pyplot.cm.jet
    n = data.shape[0] / 2

    fig1 = pyplot.figure()
    ax = fig1.add_subplot(1,1,1)
    ax.hexbin(data[:,0], data[:,1], bins=50, cmap=cmap, marginals=True)
    for mean, covar in zip(clf.means, clf.covars):
        ellipse = el.Ellipse(el.Quadrupole(covar[0,0], covar[1,1], covar[0,1]),
                             geom.Point2D(mean[0], mean[1]))
        ellipse.plot(fill=False, axes=ax, show=False)
    fig1.canvas.draw()

    import mpl_toolkits.mplot3d
    h, xe, ye = numpy.histogram2d(data[:n,0], data[:n,1], bins=(25,25), normed=True)
    xc = 0.5*(xe[1:] + xe[:-1])
    yc = 0.5*(ye[1:] + ye[:-1])
    xg, yg = numpy.meshgrid(xc, yc)

    fig2 = pyplot.figure(figsize=(18, 6))
    ax = fig2.add_subplot(1,3,1, projection='3d')
    ax.plot_surface(xg, yg, h.transpose(), rstride=1, cstride=1, shade=False, cmap=cmap,
                    linewidth=0, antialiased=False, vmin=0, vmax=0.8)
    ax.set_xlabel("ellipticity")
    ax.set_ylabel("radius")
    ax.set_zlim(0, 0.8)

    obs1 = numpy.array([xg.flatten(), yg.flatten()]).transpose()
    obs2 = numpy.array([-xg.flatten(), yg.flatten()]).transpose()
    logprob1, posteriors1 = clf.eval(obs1)
    logprob2, posteriors2 = clf.eval(obs2)
    z = (numpy.exp(logprob1) + numpy.exp(logprob2)).reshape(xg.shape)

    ax = fig2.add_subplot(1,3,2, projection='3d')
    ax.plot_surface(xg, yg, z, rstride=1, cstride=1, shade=False, cmap=cmap,
                    linewidth=0, antialiased=False, vmin=0, vmax=0.8)
    ax.set_xlabel("ellipticity")
    ax.set_ylabel("radius")
    ax.set_zlim(0, 0.8)

    ax = fig2.add_subplot(1,3,3, projection='3d')
    ax.plot_surface(xg, yg, h.transpose() - z,
                    rstride=1, cstride=1, shade=False, cmap=cmap,
                    linewidth=0, antialiased=False)
    ax.set_xlabel("ellipticity")
    ax.set_ylabel("radius")
    ax.set_zlim(-0.4, 0.4)

def transform(ellipticity, radius):
    n = ellipticity.size
    data = numpy.zeros((2*n, 2), dtype=float)
    data[:n,0] = ellipticity
    data[:n,1] = radius
    data[n:,0] = -ellipticity
    data[n:,1] = radius
    return data

if __name__ == "__main__":
    if sys.argv[1] == "convert":
        convert(sys.argv[2])
    elif sys.argv[1] == "plot":
        plot(sys.argv[2])
    elif sys.argv[1] == "fit":
        fit(sys.argv[2])

