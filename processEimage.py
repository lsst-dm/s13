import lsst.meas.multifit.starSelector

root.calibrate.astrometry.forceKnownWcs = True

# The PSF should be super-well constrained, and we know where all the stars
# are, so let's not reject any of them.
# We have enough signal to bump these numbers up even higher, but then
# things slow down considerably.
root.calibrate.measurePsf.starSelector.name = "s13"
root.calibrate.measurePsf.psfDeterminer["pca"].nEigenComponents = 8
root.calibrate.measurePsf.psfDeterminer["pca"].spatialOrder = 2
root.calibrate.measurePsf.psfDeterminer["pca"].sizeCellX = 200
root.calibrate.measurePsf.psfDeterminer["pca"].sizeCellY = 200
root.calibrate.measurePsf.psfDeterminer["pca"].nStarPerCell = 8
root.calibrate.measurePsf.psfDeterminer["pca"].nStarPerCellSpatialFit = 16
root.calibrate.measurePsf.psfDeterminer["pca"].spatialReject = 5
root.calibrate.measurePsf.psfDeterminer["pca"].nIterForPsf = 0

root.doAddNoise = True
root.noiseValue = 1000
root.measurement.algorithms["flux.sinc"].radius2 = 5.0

root.calibrate.background.algorithm = "CONSTANT"
root.calibrate.detection.background.algorithm = "CONSTANT"
root.detection.background.algorithm = "CONSTANT"
root.calibrate.repair.doCosmicRay = False
root.calibrate.repair.doInterpolate = False
