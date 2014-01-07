import lsst.meas.multifit
root.fitter.retarget(lsst.meas.multifit.AdaptiveImportanceSamplerTask)
root.fitter.doMarginalizeAmplitudes = False
root.tag = "direct"
root.previous = "marginal"
