import os

import lsst.meas.multifit
import lsst.afw.geom
import lsst.afw.table
import lsst.afw.display.ds9
import lsst.daf.persistence

class Fitter(object):

    def __init__(self, rerun="disk-0100", dataId=None, config=None):
        if dataId is None:
            dataId = dict(visit=1, raft="2,2", sensor="1,1")
        root = os.environ["S13_DATA_DIR"]
        self.butler = lsst.daf.persistence.Butler(os.path.join(root, "output", rerun))
        self.dataRef = self.butler.dataRef("calexp", dataId=dataId)
        self.exposure = self.dataRef.get("calexp", immediate=True)
        self.modelfits = self.dataRef.get("modelfits", immediate=True)
        self.task = lsst.meas.multifit.MeasureCcdTask(config=config)

    def display(self):
        lsst.afw.display.ds9.mtv(self.exposure)
        with lsst.afw.display.ds9.Buffering():
            for record in self.modelfits:
                lsst.afw.display.ds9.dot(
                    record.getShape(),
                    record.getX(),
                    record.getY()
                )
                lsst.afw.display.ds9.dot(
                    str(record.getId()),
                    record.getX(),
                    record.getY()
                )

