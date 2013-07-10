root.measurement.algorithms["flux.sinc"].radius2 = 5.0

root.calibrate.doBackground = False
root.calibrate.detection.reEstimateBackground = False
root.detection.reEstimateBackground = False
root.calibrate.repair.doCosmicRay = False
root.calibrate.repair.doInterpolate = False
root.doDeblend = False

import lsst.meas.extensions.multiShapelet
root.measurement.algorithms.names -= lsst.meas.extensions.multiShapelet.algorithms
root.measurement.slots.modelFlux = root.measurement.slots.instFlux
