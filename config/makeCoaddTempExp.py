from lsst.pipe.tasks.selectImages import WcsSelectImagesTask
root.bgSubtracted = True
root.doApplyUberCal = True
root.select.retarget(WcsSelectImagesTask)
