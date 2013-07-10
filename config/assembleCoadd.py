from lsst.pipe.tasks.selectImages import WcsSelectImagesTask
from lsst.pipe.tasks.scaleZeroPoint import ScaleZeroPointTask
root.select.retarget(WcsSelectImagesTask)
root.scaleZeroPoint.retarget(ScaleZeroPointTask)
root.doSigmaClip = False
root.doMatchBackgrounds = False
