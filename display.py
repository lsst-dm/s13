import os

import lsst.meas.multifit
import lsst.afw.display.ds9
import lsst.daf.persistence

class Fitter(object):

    def __init__(self, rerun="disk-01", dataId=None, config=None):
        if dataId is None:
            dataId = dict(visit=1, raft="2,2", sensor="1,1")
        root = os.environ["S13_DATA_DIR"]
        self.butler = lsst.daf.persistence.Butler(os.path.join(root, "output", rerun))
        self.exposure = self.butler.get("calexp", dataId, immediate=True)
        self.sources = self.butler.get("src", dataId, immediate=True)
        self.task = lsst.meas.multifit.MeasureCcdTask(config=config)
        self.output = lsst.meas.multifit.ModelFitCatalog(self.task.schema)

    def __call__(self, srcId):
        record = self.output.addNew()
        self.task.processObject(self.exposure, self.sources.find(srcId), record)
        return record

    def display(self):
        lsst.afw.display.ds9.mtv(self.exposure)
        with lsst.afw.display.ds9.Buffering():
            for source in self.sources:
                lsst.afw.display.ds9.dot(
                    source.getShape(),
                    source.getX(),
                    source.getY()
                )
                lsst.afw.display.ds9.dot(
                    str(source.getId()),
                    source.getX(),
                    source.getY()
                )

