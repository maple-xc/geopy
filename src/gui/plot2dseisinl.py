#############################################################################################
#                                                                                           #
# Author:   Haibin Di                                                                       #
# Date:     March 2018                                                                      #
#                                                                                           #
#############################################################################################

# Create a window for plot seismic inline slices


from PyQt5 import QtCore, QtGui, QtWidgets
import numpy as np
import sys, os
#
sys.path.append(os.path.dirname(__file__)[:-4])
from basic.data import data as basic_data
from core.settings import settings as core_set
from seismic.analysis import analysis as seis_ays
from seismic.visualization import visualization as seis_vis
from vis.colormap import colormap as vis_cmap
from gui.configplayer import configplayer as gui_configplayer

QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)


class plot2dseisinl(object):

    survinfo = {}
    seisdata = {}
    plotstyle = core_set.Visual['Image']
    playerconfig = core_set.Visual['Player']
    fontstyle = core_set.Visual['Font']
    #
    iconpath = os.path.dirname(__file__)
    dialog = None

    def setupGUI(self, Plot2DSeisInl):
        Plot2DSeisInl.setObjectName("Plot2DSeisInl")
        Plot2DSeisInl.setFixedSize(400, 450)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(os.path.join(self.iconpath, "icons/visinl.png")),
                       QtGui.QIcon.Normal,
                       QtGui.QIcon.Off)
        Plot2DSeisInl.setWindowIcon(icon)
        #
        self.lblattrib = QtWidgets.QLabel(Plot2DSeisInl)
        self.lblattrib.setObjectName("lblattrib")
        self.lblattrib.setGeometry(QtCore.QRect(10, 10, 150, 30))
        self.lwgattrib = QtWidgets.QListWidget(Plot2DSeisInl)
        self.lwgattrib.setObjectName("lwgattrib")
        self.lwgattrib.setGeometry(QtCore.QRect(160, 10, 230, 200))
        self.lwgattrib.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.lblslice = QtWidgets.QLabel(Plot2DSeisInl)
        self.lblslice.setObjectName("lblslice")
        self.lblslice.setGeometry(QtCore.QRect(10, 230, 150, 30))
        self.slbslice = QtWidgets.QScrollBar(Plot2DSeisInl)
        self.slbslice.setObjectName("slbslice")
        self.slbslice.setOrientation(QtCore.Qt.Horizontal)
        self.slbslice.setGeometry(QtCore.QRect(160, 230, 170, 30))
        self.ldtslice = QtWidgets.QLineEdit(Plot2DSeisInl)
        self.ldtslice.setObjectName("ldtslice")
        self.ldtslice.setGeometry(QtCore.QRect(340, 230, 50, 30))
        self.ldtslice.setAlignment(QtCore.Qt.AlignCenter)
        self.lblcmap = QtWidgets.QLabel(Plot2DSeisInl)
        self.lblcmap.setObjectName("lblcmap")
        self.lblcmap.setGeometry(QtCore.QRect(10, 270, 150, 30))
        self.cbbcmap = QtWidgets.QComboBox(Plot2DSeisInl)
        self.cbbcmap.setObjectName("cbbcmap")
        self.cbbcmap.setGeometry(QtCore.QRect(160, 270, 170, 30))
        self.cbxflip = QtWidgets.QCheckBox(Plot2DSeisInl)
        self.cbxflip.setObjectName("cbxflip")
        self.cbxflip.setGeometry(QtCore.QRect(340, 270, 50, 30))
        self.lblrange = QtWidgets.QLabel(Plot2DSeisInl)
        self.lblrange.setObjectName("lblrange")
        self.lblrange.setGeometry(QtCore.QRect(10, 310, 150, 30))
        self.ldtmin = QtWidgets.QLineEdit(Plot2DSeisInl)
        self.ldtmin.setObjectName("ldtmin")
        self.ldtmin.setGeometry(QtCore.QRect(160, 310, 90, 30))
        self.ldtmin.setAlignment(QtCore.Qt.AlignCenter)
        self.ldtmax = QtWidgets.QLineEdit(Plot2DSeisInl)
        self.ldtmax.setObjectName("ldtmax")
        self.ldtmax.setGeometry(QtCore.QRect(300, 310, 90, 30))
        self.ldtmax.setAlignment(QtCore.Qt.AlignCenter)
        self.lblrangeto = QtWidgets.QLabel(Plot2DSeisInl)
        self.lblrangeto.setObjectName("lblrangeto")
        self.lblrangeto.setGeometry(QtCore.QRect(250, 310, 50, 30))
        self.lblrangeto.setAlignment(QtCore.Qt.AlignCenter)
        #
        self.btnconfigplayer = QtWidgets.QPushButton(Plot2DSeisInl)
        self.btnconfigplayer.setObjectName("btnconfigplayer")
        self.btnconfigplayer.setGeometry(QtCore.QRect(230, 350, 160, 30))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(os.path.join(self.iconpath, "icons/video.png")),
                       QtGui.QIcon.Normal,
                       QtGui.QIcon.Off)
        self.btnconfigplayer.setIcon(icon)
        #
        self.btnplotslice = QtWidgets.QPushButton(Plot2DSeisInl)
        self.btnplotslice.setObjectName("btnplotslice")
        self.btnplotslice.setGeometry(QtCore.QRect(120, 400, 160, 30))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(os.path.join(self.iconpath, "icons/visinl.png")),
                       QtGui.QIcon.Normal,
                       QtGui.QIcon.Off)
        self.btnplotslice.setIcon(icon)
        #
        self.msgbox = QtWidgets.QMessageBox(Plot2DSeisInl)
        self.msgbox.setObjectName("msgbox")
        _center_x = Plot2DSeisInl.geometry().center().x()
        _center_y = Plot2DSeisInl.geometry().center().y()
        self.msgbox.setGeometry(QtCore.QRect(_center_x - 150, _center_y - 50, 300, 100))
        #
        self.retranslateGUI(Plot2DSeisInl)
        QtCore.QMetaObject.connectSlotsByName(Plot2DSeisInl)



    def retranslateGUI(self, Plot2DSeisInl):
        self.dialog = Plot2DSeisInl
        #
        _translate = QtCore.QCoreApplication.translate
        Plot2DSeisInl.setWindowTitle(_translate("Plot2DSeisInl", "2D Window: Seismic Inline"))
        self.lblattrib.setText(_translate("Plot2DSeisInl", "Select target properties:"))
        if (self.checkSurvInfo() is True) and (self.checkSeisData() is True):
            _firstattrib = None
            for i in sorted(self.seisdata.keys()):
                item = QtWidgets.QListWidgetItem(self.lwgattrib)
                item.setText(_translate("Plot2DSeisInl", i))
                self.lwgattrib.addItem(item)
                if _firstattrib is None:
                    _firstattrib = item
            self.lwgattrib.setCurrentItem(_firstattrib)
        self.lwgattrib.itemSelectionChanged.connect(self.changeLwgAttrib)
        #
        self.lblslice.setText(_translate("Plot2DSeisInl", "Select inline slice:"))
        if (self.checkSurvInfo() is True):
            _slices = self.survinfo['ILRange'].astype(int)
            _slicemin = np.min(_slices)
            _slicemax = np.max(_slices)
        else:
            _slicemin = 0
            _slicemax = 0
        self.slbslice.setMinimum(_slicemin)
        self.slbslice.setMaximum(_slicemax)
        self.slbslice.setValue(_slicemin)
        if (self.checkSurvInfo() is True):
            self.ldtslice.setText(_translate("Plot2DSeisInl", str(_slicemin)))
        else:
            self.ldtslice.setText(_translate("Plot2DSeisInl", ''))
        self.slbslice.valueChanged.connect(self.changeSlbSlice)
        self.ldtslice.textChanged.connect(self.changeLdtSlice)
        #
        self.lblcmap.setText(_translate("Plot2DXlSlice", "\t  Color map:"))
        self.cbbcmap.addItems(vis_cmap.ColorMapList)
        for _i in range(len(vis_cmap.ColorMapList)):
            self.cbbcmap.setItemIcon(_i, QtGui.QIcon(
                QtGui.QPixmap(os.path.join(self.iconpath, "icons/cmap_" + vis_cmap.ColorMapList[_i] + ".png")).scaled(80, 30)))
        self.cbbcmap.setCurrentIndex(list.index(vis_cmap.ColorMapList, self.plotstyle['Colormap']))
        #
        self.cbxflip.setText(_translate("Plot2DSeisInl", ""))
        self.cbxflip.setIcon(QtGui.QIcon(
            QtGui.QPixmap(os.path.join(self.iconpath, "icons/flip.png")).scaled(80, 80)))
        #
        self.lblrange.setText(_translate("Plot2DSeisInl", "\t       Range:"))
        self.lblrangeto.setText(_translate("Plot2DSeisInl", "~~~"))
        if (self.checkSurvInfo() is True) \
                and (self.checkSeisData() is True) \
                and (self.lwgattrib.currentItem() is not None):
            _min, _max = self.getAttribRange(self.lwgattrib.currentItem().text())
            self.ldtmin.setText(_translate("Plot2DSeisInl", str(_min)))
            self.ldtmax.setText(_translate("Plot2DSeisInl", str(_max)))
        #
        self.btnconfigplayer.setText(_translate("Plot2DSeisInl", "Player Configuration"))
        self.btnconfigplayer.clicked.connect(self.clickBtnConfigPlayer)
        #
        self.btnplotslice.setText(_translate("Plot2DSeisInl", "Seismic Inline Viewer"))
        self.btnplotslice.setDefault(True)
        self.btnplotslice.clicked.connect(self.clickBtnPlotSlice)


    def changeLwgAttrib(self):
        if len(self.lwgattrib.selectedItems()) > 0:
            # if len(self.lwgattrib.selectedItems()) > 1:
            #     self.ldtmin.setEnabled(False)
            #     self.ldtmax.setEnabled(False)
            # else:
            #     self.ldtmin.setEnabled(True)
            #     self.ldtmax.setEnabled(True)
            _slices = self.survinfo['ILRange'].astype(int)
            _slicemin = np.min(_slices)
            _slicemax = np.max(_slices)
            self.slbslice.setMinimum(_slicemin)
            self.slbslice.setMaximum(_slicemax)
            self.ldtslice.setText(str(self.slbslice.value()))
            _min, _max = self.getAttribRange(self.lwgattrib.currentItem().text())
            self.ldtmin.setText(str(_min))
            self.ldtmax.setText(str(_max))
        else:
            self.slbslice.setMinimum(0)
            self.slbslice.setMaximum(0)
            self.ldtslice.setText('')
            self.ldtmin.setText('')
            self.ldtmax.setText('')


    def changeSlbSlice(self):
        self.ldtslice.setText(str(self.slbslice.value()))


    def changeLdtSlice(self):
        if len(self.ldtslice.text()) > 0:
            _val = basic_data.str2int(self.ldtslice.text())
            if _val >= self.slbslice.minimum() and _val <= self.slbslice.maximum():
                self.slbslice.setValue(_val)


    def clickBtnConfigPlayer(self):
        _config = QtWidgets.QDialog()
        _gui = gui_configplayer()
        _gui.playerconfig = {}
        _gui.playerconfig['Inline'] = self.playerconfig
        _gui.setupGUI(_config)
        _config.exec()
        self.playerconfig = _gui.playerconfig['Inline']
        _config.show()


    def clickBtnPlotSlice(self):
        self.refreshMsgBox()
        #
        _attriblist = self.lwgattrib.selectedItems()
        if len(_attriblist) < 1:
            print("Plot2DSeisInl: No property selected for plot")
            QtWidgets.QMessageBox.critical(self.msgbox,
                                           '2D Window: Seismic Inline',
                                           'No property selected for plot')
            return
        #
        _sls = self.slbslice.value()
        _cmap = self.cbbcmap.currentIndex()
        _flip = self.cbxflip.isChecked()
        _min = basic_data.str2float(self.ldtmin.text())
        _max = basic_data.str2float(self.ldtmax.text())
        if _min is False or _max is False:
            print("Plot2DSeisInl: Non-float range specified for plot")
            QtWidgets.QMessageBox.critical(self.msgbox,
                                           '2D Window: Seismic Inline',
                                           'Non-float range specified for plot')
            return
        #
        for i in range(len(_attriblist)):
            print("Plot2DSeisInl: Plot %d of %d properties: %s" % (i + 1, len(_attriblist), _attriblist[i].text()))
            _data = self.seisdata[_attriblist[i].text()]
            _data = np.reshape(_data, [_data.shape[0], -1])
            _data = np.mean(_data, axis=1)
            _data = np.reshape(_data, [len(_data), 1])
            _data = np.transpose(np.reshape(_data, [self.survinfo['ILNum'], self.survinfo['XLNum'], self.survinfo['ZNum']]),
                                 [2, 1, 0])
            # if len(_attriblist) > 1:
            #     _min, _max = self.getAttribRange(_attriblist[i].text())
            seis_vis.plotSeisILSlicePlayerFrom3DMat(_data, initinlsl=_sls, seisinfo=self.survinfo,
                                                    titlesurf=': ' + _attriblist[i].text(),
                                                    valuemin=_min,
                                                    valuemax=_max,
                                                    colormap=vis_cmap.ColorMapList[_cmap],
                                                    flipcmap=_flip, colorbaron=True,
                                                    interpolation=self.plotstyle['Interpolation'].lower(),
                                                    playerconfig=self.playerconfig,
                                                    fontstyle=self.fontstyle,
                                                    qicon=QtGui.QIcon(os.path.join(self.iconpath, "icons/logo.png"))
                                                    )
        #
        # QtWidgets.QMessageBox.information(self.msgbox,
        #                                   "Plot Inline",
        #                                   str(len(_attriblist)) + " properties plotted successfully")
        # self.dialog.close()
        return


    def getAttribRange(self, f):
        _min = -1
        _max = 1
        if (self.checkSurvInfo() is True) \
                and (self.checkSeisData() is True) \
                and (f in self.seisdata.keys()):
            _min = np.min(self.seisdata[f])
            _max = np.max(self.seisdata[f])
        return _min, _max


    def refreshMsgBox(self):
        _center_x = self.dialog.geometry().center().x()
        _center_y = self.dialog.geometry().center().y()
        self.msgbox.setGeometry(QtCore.QRect(_center_x - 150, _center_y - 50, 300, 100))


    def checkSurvInfo(self):
        self.refreshMsgBox()
        #
        if seis_ays.checkSeisInfo(self.survinfo) is False:
            # print("Plot2DSeisInl: Survey not found")
            # QtWidgets.QMessageBox.critical(self.msgbox,
            #                                'Plot Seismic Inline',
            #                                'Survey not found')
            return False
        return True


    def checkSeisData(self):
        self.refreshMsgBox()
        #
        for f in self.seisdata.keys():
            if np.shape(self.seisdata[f])[0] != self.survinfo['SampleNum']:
                # print("Plot2DSeisInl: Seismic & survey not match")
                # QtWidgets.QMessageBox.critical(self.msgbox,
                #                                'Plot Seismic Inline',
                #                                'Seismic & survey not match')
                return False
        return True


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    Plot2DSeisInl = QtWidgets.QWidget()
    gui = plot2dseisinl()
    gui.setupGUI(Plot2DSeisInl)
    Plot2DSeisInl.show()
    sys.exit(app.exec_())