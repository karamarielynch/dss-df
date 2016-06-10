import sys
from PyQt4 import QtGui, QtCore
import pyqtgraph as pg
import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import odrpack as odr
import scipy.integrate as integrate
import lineshape
from datetime import datetime


class FitDSSData(QtGui.QMainWindow):

    def __init__(self):
        super(FitDSSData, self).__init__()

        self.initUI()


    def initUI(self):

        self.centralWidget = QtGui.QWidget()
        self.setCentralWidget(self.centralWidget)
        self.grid = QtGui.QGridLayout(self.centralWidget)

        gridLayout = QtGui.QGridLayout()
        self.grid.addLayout(gridLayout,0,0,1,1)
        self.grid.addWidget(QtGui.QWidget(),1,0)

        self.setWindowTitle('DSS Data Fitter')

        self.tabs = QtGui.QTabWidget()

        # PLOT TAB

        self.tab_plot = QtGui.QWidget()
        plot_grid = QtGui.QGridLayout()

        configFile = sys.argv[1]
        self.configData = np.genfromtxt(configFile, delimiter='\t', dtype="S8,S8,S8,S8,f8,f8", names=['Det', 'Name', 'Type', 'Number', 'Slope', 'Offset'] )
        self.Dets, self.Names, self.Types, self.Numbers, self.Slopes, self.Offsets  = self.configData["Det"].astype(int) - 140, self.configData["Name"].astype(str), self.configData["Type"].astype(int), self.configData["Number"].astype(int), self.configData["Slope"], self.configData["Offset"]
        self.Dets, self.Names, self.Types, self.Numbers, self.Slopes, self.Offsets  = self.Dets[(self.Names != "empty")], self.Names[(self.Names != "empty")], self.Types[(self.Names != "empty")], self.Numbers[(self.Names != "empty")], self.Slopes[(self.Names != "empty")], self.Offsets[(self.Names != "empty")]

        DetType = QtGui.QLabel('Detector Type')
        plot_grid.addWidget(DetType, 0, 0)
        DetTypeEdit = QtGui.QLabel()
        self.DetTypeList = QtGui.QComboBox(self)
        self.DetTypeList.addItems(self.Names)
        plot_grid.addWidget(self.DetTypeList, 0, 1)
        self.DetTypeList.activated[str].connect(self.updateIndex)
        self.DetTypeList.activated[str].connect(self.calibrateData)
        self.DetTypeList.activated[str].connect(self.reformatData)
        self.DetTypeList.activated[str].connect(self.updatePlotGate)
        self.DetTypeList.activated[str].connect(self.updatePlot)
        self.DetTypeList.activated[str].connect(self.calcRate)

        LoadButton = QtGui.QPushButton('Load Data', self)
        plot_grid.addWidget(LoadButton, 0, 3)
        LoadButton.clicked.connect(self.chooseFile)
        LoadButton.clicked.connect(self.loadData)
        LoadButton.clicked.connect(self.updateIndex)
        LoadButton.clicked.connect(self.calibrateData)
        LoadButton.clicked.connect(self.reformatData)
        LoadButton.clicked.connect(self.updatePlot)
        LoadButton.clicked.connect(self.calcRate)
        LoadButton.clicked.connect(self.findCoincs)
        LoadButton.setToolTip('Load data')
        self.loadDataBoolean = False

        Parameter   = QtGui.QLabel('Plot')
        plot_grid.addWidget(Parameter, 1, 0)
        ParameterEdit = QtGui.QLabel()
        self.ParameterList = QtGui.QComboBox(self)
        self.ParameterNames = ["Energy", "Channel", "Timestamp"]
        self.ParameterList.addItems(self.ParameterNames)
        plot_grid.addWidget(self.ParameterList, 1, 1)
        self.ParameterList.activated[str].connect(self.updatePlotGate)
        self.ParameterList.activated[str].connect(self.reformatData)
        self.ParameterList.activated[str].connect(self.updatePlot)

        NoBins  = QtGui.QLabel('Bins')
        plot_grid.addWidget(NoBins, 1, 2)
        self.NoBinsEdit = pg.SpinBox(value=1500, dec=True, minStep=1)
        self.NoBinsEdit.setRange(1, 10**8)
        plot_grid.addWidget(self.NoBinsEdit, 1, 3)
        self.NoBinsEdit.valueChanged.connect(self.reformatData)
        self.NoBinsEdit.valueChanged.connect(self.updatePlot)

        PlotFrom    = QtGui.QLabel('Plot: from')
        plot_grid.addWidget(PlotFrom, 2, 0)
        self.PlotFromEdit = pg.SpinBox(value=0, dec=True, minStep=1)
        plot_grid.addWidget(self.PlotFromEdit, 2, 1)
        self.PlotFromEdit.valueChanged.connect(self.reformatData)
        self.PlotFromEdit.valueChanged.connect(self.updatePlot)

        PlotTo= QtGui.QLabel('to')
        plot_grid.addWidget(PlotTo, 2, 2)
        self.PlotToEdit = pg.SpinBox(value=1500, dec=True, minStep=1)
        plot_grid.addWidget(self.PlotToEdit, 2, 3)
        self.PlotToEdit.valueChanged.connect(self.reformatData)
        self.PlotToEdit.valueChanged.connect(self.updatePlot)

        Calibration_m   = QtGui.QLabel('Calibration: slope')
        plot_grid.addWidget(Calibration_m, 3, 0)
        self.Calibration_mEdit = pg.SpinBox(value=self.Slopes[0])
        self.Calibration_mEdit.setRange(0, 100)
        plot_grid.addWidget(self.Calibration_mEdit, 3, 1)
        self.Calibration_mEdit.valueChanged.connect(self.reformatData)
        self.Calibration_mEdit.valueChanged.connect(self.updatePlot)

        Calibration_c   = QtGui.QLabel('offset')
        plot_grid.addWidget(Calibration_c, 3, 2)
        self.Calibration_cEdit = pg.SpinBox(value=self.Offsets[0])
        self.Calibration_cEdit.setRange(-100, 100)
        plot_grid.addWidget(self.Calibration_cEdit, 3, 3)
        self.Calibration_cEdit.valueChanged.connect(self.reformatData)
        self.Calibration_cEdit.valueChanged.connect(self.updatePlot)


        self.tab_plot.setLayout(plot_grid)

        # FIT TAB

        self.tab_fit = QtGui.QWidget()
        fit_grid = QtGui.QGridLayout()

        LineProfile = QtGui.QLabel('Line Profile')
        fit_grid.addWidget(LineProfile, 0, 0)
        LineProfileEdit = QtGui.QLabel()
        self.LineProfileList = QtGui.QComboBox(self)
        self.LineProfileNames = ["Gaussian", "CrystalBall"]
        self.LineProfileList.addItems(self.LineProfileNames)
        fit_grid.addWidget(self.LineProfileList, 0, 1)
        self.LineProfileList.activated[str].connect(self.updatePlot)

        NoPeaks = QtGui.QLabel('Peaks')
        fit_grid.addWidget(NoPeaks, 0, 2)
        NoPeaksEdit = QtGui.QLabel()
        self.NoPeaksList = QtGui.QComboBox(self)
        self.NoPeaksNumbers = ["1", "2"]
        self.NoPeaksList.addItems(self.NoPeaksNumbers)
        fit_grid.addWidget(self.NoPeaksList, 0, 3)
        self.NoPeaksList.activated[str].connect(self.updateNoPeaks)
        self.NoPeaksList.activated[str].connect(self.updatePlot)

        self.PeaksButton = QtGui.QPushButton('Show Peaks', self)
        fit_grid.addWidget(self.PeaksButton, 0, 4)
        self.PeaksButton.clicked.connect(self.updatePlot)
        self.PeaksButton.setToolTip('Plot individual peaks')
        self.PeaksButton.setCheckable(True)

        FitButton = QtGui.QPushButton('Fit Data', self)
        FitButton.clicked.connect(self.fitData)
        FitButton.setToolTip('Fit data')
        fit_grid.addWidget(FitButton, 0, 5)

        start   = QtGui.QLabel('Start (keV)')
        fit_grid.addWidget(start, 1, 0)
        self.startEdit = pg.SpinBox(value=5000, dec=True, minStep=1)
        self.startEdit.setRange(0, 10000)
        fit_grid.addWidget(self.startEdit, 1, 1)
        self.startEdit.valueChanged.connect(self.updatePlot)

        stop    = QtGui.QLabel('Stop (keV)')
        fit_grid.addWidget(stop, 1, 2)
        self.stopEdit = pg.SpinBox(value=5250, dec=True, minStep=1)
        self.stopEdit.setRange(0, 10000)
        fit_grid.addWidget(self.stopEdit, 1, 3)
        self.stopEdit.valueChanged.connect(self.updatePlot)

        peak1       = QtGui.QLabel('Peak 1')
        fit_grid.addWidget(peak1, 3, 0)

        peak2       = QtGui.QLabel('Peak 2')
        fit_grid.addWidget(peak2, 4, 0)

        peak    = QtGui.QLabel('Peak (keV)')
        fit_grid.addWidget(peak, 2, 1)

        intensity   = QtGui.QLabel('Intensity')
        fit_grid.addWidget(intensity, 2, 2)

        FWHM    = QtGui.QLabel('FWHM (keV)')
        fit_grid.addWidget(FWHM, 2, 3)

        background  = QtGui.QLabel('Background')
        fit_grid.addWidget(background, 2, 4)

        content = QtGui.QLabel('Integral value')
        fit_grid.addWidget(content, 2, 5)

        redchi2 = QtGui.QLabel('Red. chi2')
        fit_grid.addWidget(redchi2, 5, 0)

        peak1_params = np.array([5105, 200, 25, 0])
        self.peak1_paramsEdits = []
        for i in range(4):
            self.peak1_paramsEdit = pg.SpinBox(value=0, dec=True, minStep=1)
            self.peak1_paramsEdit.setRange(0, 10**6)
            fit_grid.addWidget(self.peak1_paramsEdit, 3, 1+i)
            self.peak1_paramsEdits.append(self.peak1_paramsEdit)
            self.peak1_paramsEdit.setValue(peak1_params[i])
            self.peak1_paramsEdit.valueChanged.connect(self.updatePlot)

        peak2_params = np.array([5155, 2000, 30])
        self.peak2_paramsEdits = []
        for i in range(3):
            self.peak2_paramsEdit = pg.SpinBox(value=0, dec=True, minStep=1)
            self.peak2_paramsEdit.setRange(0, 10**6)
            fit_grid.addWidget(self.peak2_paramsEdit, 4, 1+i)
            self.peak2_paramsEdits.append(self.peak2_paramsEdit)
            self.peak2_paramsEdit.setValue(peak2_params[i])
            self.peak2_paramsEdit.valueChanged.connect(self.updatePlot)
            self.peak2_paramsEdit.setEnabled(False)

        self.peak1_contentEdit = QtGui.QLabel('0')
        fit_grid.addWidget(self.peak1_contentEdit, 3, 5)

        self.peak2_contentEdit = QtGui.QLabel('0')
        fit_grid.addWidget(self.peak2_contentEdit, 4, 5)

        self.peak1_chi2Edit = QtGui.QLabel('0')
        fit_grid.addWidget(self.peak1_chi2Edit, 5, 1)

        self.tab_fit.setLayout(fit_grid)


        # RATE TAB

        self.tab_rate = QtGui.QWidget()
        rate_grid = QtGui.QGridLayout()

        TotalTime   = QtGui.QLabel('Total time')
        rate_grid.addWidget(TotalTime, 0, 0)
        self.TotalTimeEdit = QtGui.QLabel('0')
        rate_grid.addWidget(self.TotalTimeEdit, 0, 1)

        TotalCounts = QtGui.QLabel('Total counts')
        rate_grid.addWidget(TotalCounts, 0, 2)
        self.TotalCountsEdit = QtGui.QLabel('0')
        rate_grid.addWidget(self.TotalCountsEdit, 0, 3)

        CountRate   = QtGui.QLabel('Count rate')
        rate_grid.addWidget(CountRate, 0, 4)
        self.CountRateEdit = QtGui.QLabel('0')
        rate_grid.addWidget(self.CountRateEdit, 0, 5)

        Gate_i  = QtGui.QLabel('Energy gate: from')
        rate_grid.addWidget(Gate_i, 1, 0)
        self.Gate_iEdit = pg.SpinBox(value=0, dec=True, minStep=1)
        rate_grid.addWidget(self.Gate_iEdit, 1, 1)
        self.Gate_iEdit.valueChanged.connect(self.calcGatedRate)

        Gate_f  = QtGui.QLabel('to')
        rate_grid.addWidget(Gate_f, 1, 2)
        self.Gate_fEdit = pg.SpinBox(value=0, dec=True, minStep=1)
        rate_grid.addWidget(self.Gate_fEdit, 1, 3)
        self.Gate_fEdit.valueChanged.connect(self.calcGatedRate)

        GatedCounts = QtGui.QLabel('Gated counts')
        rate_grid.addWidget(GatedCounts, 1, 4)
        self.GatedCountsEdit = QtGui.QLabel('0')
        rate_grid.addWidget(self.GatedCountsEdit, 1, 5)

        GatedRate   = QtGui.QLabel('Gated rate')
        rate_grid.addWidget(GatedRate, 1, 6)
        self.GatedRateEdit = QtGui.QLabel('0')
        rate_grid.addWidget(self.GatedRateEdit, 1, 7)

        self.tab_rate.setLayout(rate_grid)


        # COINC TAB

        self.tab_coinc = QtGui.QWidget()
        coinc_grid = QtGui.QGridLayout()

        TimeStart   = QtGui.QLabel('Time window: Start')
        coinc_grid.addWidget(TimeStart, 0, 0)
        self.TimeStartEdit = pg.SpinBox(value=0)
        coinc_grid.addWidget(self.TimeStartEdit, 0, 1)
        self.TimeStartEdit.valueChanged.connect(self.recalcCoincs)

        TimeStop    = QtGui.QLabel('Stop')
        coinc_grid.addWidget(TimeStop, 0, 2)
        self.TimeStopEdit = pg.SpinBox(value=15)
        coinc_grid.addWidget(self.TimeStopEdit, 0, 3)
        self.TimeStopEdit.valueChanged.connect(self.recalcCoincs)

        DetectorPlot    = QtGui.QLabel('Plot detector')
        coinc_grid.addWidget(DetectorPlot, 1, 0)
        self.DetectorPlotList = QtGui.QComboBox(self)
        self.DetectorPlotList.addItems(self.Names)
        coinc_grid.addWidget(self.DetectorPlotList, 1, 1)
        self.DetectorPlotList.activated[str].connect(self.recalcCoincs)

        NoBinsPlot  = QtGui.QLabel('Bins')
        coinc_grid.addWidget(NoBinsPlot, 1, 2)
        self.NoBinsPlotEdit = pg.SpinBox(value=100, dec=True, minStep=1)
        self.NoBinsPlotEdit.setRange(1, 10**8)
        coinc_grid.addWidget(self.NoBinsPlotEdit, 1, 3)
        self.NoBinsPlotEdit.valueChanged.connect(self.recalcCoincs)

        DetectorGate    = QtGui.QLabel('Gate on detector')
        coinc_grid.addWidget(DetectorGate, 2, 0)
        self.DetectorGateList = QtGui.QComboBox(self)
        self.DetectorGateList.addItems(self.Names)
        coinc_grid.addWidget(self.DetectorGateList, 2, 1)
        self.DetectorGateList.setCurrentIndex(1)
        self.DetectorGateList.activated[str].connect(self.recalcCoincs)

        NoBinsGate  = QtGui.QLabel('Bins')
        coinc_grid.addWidget(NoBinsGate, 2, 2)
        self.NoBinsGateEdit = pg.SpinBox(value=100, dec=True, minStep=1)
        self.NoBinsGateEdit.setRange(1, 10**8)
        coinc_grid.addWidget(self.NoBinsGateEdit, 2, 3)

        EnergyGate_i    = QtGui.QLabel('Energy gate: from')
        coinc_grid.addWidget(EnergyGate_i, 3, 0)
        self.EnergyGate_iEdit = pg.SpinBox(value=0)
        coinc_grid.addWidget(self.EnergyGate_iEdit, 3, 1)
        self.EnergyGate_iEdit.valueChanged.connect(self.recalcCoincs)

        EnergyGate_f    = QtGui.QLabel('to')
        coinc_grid.addWidget(EnergyGate_f, 3, 2)
        self.EnergyGate_fEdit = pg.SpinBox(value=0)
        coinc_grid.addWidget(self.EnergyGate_fEdit, 3, 3)
        self.EnergyGate_fEdit.valueChanged.connect(self.recalcCoincs)

        EventsCoinc = QtGui.QLabel('Events')
        coinc_grid.addWidget(EventsCoinc, 3, 4)
        self.EventsCoincEdit = QtGui.QLabel('0')
        coinc_grid.addWidget(self.EventsCoincEdit, 3, 5)

        HistoButton = QtGui.QPushButton('Show 2D hist', self)
        coinc_grid.addWidget(HistoButton, 2, 5)
        HistoButton.clicked.connect(self.show2Dhist)
        HistoButton.setToolTip('Show 2D histogram of coincidences')

        DeltaTButton = QtGui.QPushButton("Show "+u"\u0394"+"t hist", self)
        coinc_grid.addWidget(DeltaTButton, 0, 5)
        DeltaTButton.clicked.connect(self.showdthist)
        DeltaTButton.setToolTip('Show histogram of delta t')

        self.tab_coinc.setLayout(coinc_grid)


        # ADDING WIDGETS TO TABWIDGET

        self.tabs.addTab(self.tab_plot, "Plot Data")
        self.tabs.addTab(self.tab_fit,"Fit Data")
        self.tabs.addTab(self.tab_rate,"Data Rate")
        self.tabs.addTab(self.tab_coinc,"Coincidences")
        self.grid.addWidget(self.tabs,4,0,5,11)

        self.plotWidget = pg.PlotWidget()
        self.grid.addWidget(self.plotWidget, 1, 0, 3, 11)
        self.plotWidget.setLabel('bottom', "Energy", units='eV')
        self.plotWidget.setLabel('left', "Counts per bin")
        self.estPen = pg.mkPen(color=(0, 183, 255), width = 2)

        self.vLine = pg.InfiniteLine(angle=90, movable=False, pen=(0, 183, 255))
        self.hLine = pg.InfiniteLine(angle=0, movable=False, pen=(0, 183, 255))
        self.plotWidget.addItem(self.vLine, ignoreBounds=True)
        self.plotWidget.addItem(self.hLine, ignoreBounds=True)

        self.label  = QtGui.QLabel("x = %.1f, y = %.1f" % (0, 0))
        self.grid.addWidget(self.label, 0, 10)

        self.proxy = pg.SignalProxy(self.plotWidget.scene().sigMouseMoved, rateLimit=100, slot=self.mouseMoved)

        self.show()

    def updateIndex(self):
        self.DetIndex = self.DetTypeList.currentIndex()


    def chooseFile(self):
        self.filename = QtGui.QFileDialog.getOpenFileName(self, 'Load data file')
        self.loadDataBoolean = True

    def loadData(self):
        self.startTime = datetime.now()
        np.random.seed(0)
        with open(self.filename) as f:
            length = len([line for line in f.readlines() if (len(line) > 2 and str(line[-3])=="4" and str(line[-5:-1]) != "e")])
            self.timestamps = np.ones(length) * -1
            self.channels = np.ones(length) * -1
            self.detectors = np.ones(length).astype(int)
        with open(self.filename) as f:
            for i, line in enumerate([line for line in f.readlines() if (len(line) > 2 and str(line[-3])=="4" and str(line[-5:-1]) != "e")]):
                cut_line = line[8:]
                channel = int(cut_line[:4], 16)
                timestamp = int(''.join(cut_line[5:-20].split(' ')), 16)
                detector = int(line[-2])
                if (timestamp <= 1.1*10**9): # Without this timecut, there are strange events at 8*10**12 for all data files in a very high channel [overflow?]
                    self.channels[i] = channel
                    self.timestamps[i] = timestamp
                    self.detectors[i] = detector
        self.channels, self.timestamps, self.detectors = self.channels[self.channels > -1], self.timestamps[self.channels > -1], self.detectors[self.channels > -1]
        if len(self.Slopes)>1:
            self.energies = self.channels*self.Slopes[self.detectors] + self.Offsets[self.detectors] + (np.random.rand(len(self.channels))-0.5)
        else:
            self.energies = self.channels*self.Slopes[0] + self.Offsets[0] + (np.random.rand(len(self.channels))-0.5)
        self.data = np.column_stack((self.channels, self.timestamps, self.detectors, self.energies))

    def calibrateData(self):
        m = self.Slopes[self.DetIndex]
        c = self.Offsets[self.DetIndex]
        self.Calibration_mEdit.setValue(m)
        self.Calibration_cEdit.setValue(c)

    def reformatData(self):
        detector = self.Dets[self.DetIndex]
        data    = self.data[(self.data[:,2] == detector)]
        self.channels   = data[:,0]
        self.channels   = np.array(self.channels, dtype = np.float)
        self.timestamps     = data[:,1]
        self.timestamps     = np.array(self.timestamps, dtype = np.float)


        if self.ParameterList.currentText() == "Channel":
            self.y, x = np.histogram(self.channels, bins = self.NoBinsEdit.value(), range = (float(self.PlotFromEdit.value()), float(self.PlotToEdit.value())) )
            self.x = x[0] + np.diff(x).cumsum()
            self.plotWidget.setLabel('bottom', "Channel", units='')
        if self.ParameterList.currentText() == "Timestamp":
            self.y, x = np.histogram(self.timestamps,  bins = self.NoBinsEdit.value(), range = (float(self.PlotFromEdit.value()), float(self.PlotToEdit.value())) )
            self.x = x[0] + np.diff(x).cumsum()
            self.plotWidget.setLabel('bottom', u"\u0394t", units='')
        if self.ParameterList.currentText() == "Energy":
            np.random.seed(0)
            self.channels = np.array(self.channels)
            self.energies = self.channels*self.Calibration_mEdit.value() + self.Calibration_cEdit.value() + (np.random.rand(len(self.channels))-0.5)
            self.y, x = np.histogram(self.energies,  bins = self.NoBinsEdit.value(), range = (float(self.PlotFromEdit.value()), float(self.PlotToEdit.value())) )
            self.x = x[0] + np.diff(x).cumsum()
            self.plotWidget.setLabel('bottom', "Energy", units='eV')

    def fitData(self):
        self.x_fit = np.linspace(self.startEdit.value(), self.stopEdit.value(), 1000)
        plot_mask = np.bitwise_and(self.startEdit.value() <= self.x, self.x <= self.stopEdit.value())
        self.y_cut, self.x_cut = self.y[plot_mask], self.x[plot_mask]
        if self.NoPeaksList.currentText() == "1":
            if self.LineProfileList.currentText() == "Gaussian":
                self.fitSingleGaussianPeak()
            else:
                self.fitSingleCrystalBallPeak()
        if self.NoPeaksList.currentText() == "2":
            if self.LineProfileList.currentText() == "Gaussian":
                self.fitDoubleGaussianPeak()
            else:
                self.fitDoubleCrystalBallPeak()


    def updateNoPeaks(self):
        if self.NoPeaksList.currentText() == "1":
            for i in range(3):
                self.peak2_paramsEdits[i].setEnabled(False)
        if self.NoPeaksList.currentText() == "2":
            for i in range(3):
                self.peak2_paramsEdits[i].setEnabled(True)


    def updatePlotGate(self):
        if self.ParameterList.currentText() != "Timestamp":
            if self.Types[self.DetIndex] == "0":
                self.PlotFromEdit.setValue(0.)
                self.PlotToEdit.setValue(1500.)
            if self.Types[self.DetIndex] == "1":
                self.PlotFromEdit.setValue(0.)
                self.PlotToEdit.setValue(9000.)
        else:
            self.PlotFromEdit.setValue(0.)
            self.PlotToEdit.setValue(2*10**9)


    def updatePlot(self):
        self.plotWidget.plot(clear=True)
        self.plotWidget.addItem(self.vLine, ignoreBounds=True)
        self.plotWidget.addItem(self.hLine, ignoreBounds=True)

        if  self.loadDataBoolean == True:
            histo = pg.PlotCurveItem(self.x,self.y[:-1],stepMode=True)
            self.plotWidget.addItem(histo)

        if self.ParameterList.currentText() != "Timestamp":

            self.x_fit = np.linspace(self.startEdit.value(), self.stopEdit.value(), 1000)

            if self.NoPeaksList.currentText() == "1":
                if self.LineProfileList.currentText() == "Gaussian":
                    params = np.array([self.peak1_paramsEdits[0].value(), self.peak1_paramsEdits[1].value(), self.peak1_paramsEdits[2].value(), self.peak1_paramsEdits[3].value() ])
                    self.plotWidget.plot(self.x_fit, lineshape.Single_Gaussian(params, self.x_fit), pen=self.estPen)
                else:
                    params = np.array([self.peak1_paramsEdits[0].value(), self.peak1_paramsEdits[1].value(), self.peak1_paramsEdits[2].value(), self.peak1_paramsEdits[3].value(), 0.9, 4 ])
                    self.plotWidget.plot(self.x_fit, lineshape.Single_Crystalball(params, self.x_fit), pen=self.estPen)


            if self.NoPeaksList.currentText() == "2":
                if self.PeaksButton.isChecked() == True:
                    indivPen = pg.mkPen(color=(0, 183, 255), width = 2, style=QtCore.Qt.DashLine)
                    if self.LineProfileList.currentText() == "Gaussian":
                        params1 = np.array([self.peak1_paramsEdits[0].value(), self.peak1_paramsEdits[1].value(), self.peak1_paramsEdits[2].value(), self.peak1_paramsEdits[3].value() ])
                        params2 = np.array([self.peak2_paramsEdits[0].value(), self.peak2_paramsEdits[1].value(), self.peak2_paramsEdits[2].value(), self.peak1_paramsEdits[3].value() ])
                        self.plotWidget.plot(self.x_fit, lineshape.Single_Gaussian(params1, self.x_fit), pen=indivPen)
                        self.plotWidget.plot(self.x_fit, lineshape.Single_Gaussian(params2, self.x_fit), pen=indivPen)
                    else:
                        params1 = np.array([self.peak1_paramsEdits[0].value(), self.peak1_paramsEdits[1].value(), self.peak1_paramsEdits[2].value(), self.peak1_paramsEdits[3].value(), 0.9, 4 ])
                        params2 = np.array([self.peak2_paramsEdits[0].value(), self.peak2_paramsEdits[1].value(), self.peak2_paramsEdits[2].value(), self.peak1_paramsEdits[3].value(), 0.9, 4 ])
                        self.plotWidget.plot(self.x_fit, lineshape.Single_Crystalball(params1, self.x_fit), pen=indivPen)
                        self.plotWidget.plot(self.x_fit, lineshape.Single_Crystalball(params2, self.x_fit), pen=indivPen)

                if self.PeaksButton.isChecked() == False:
                    if self.LineProfileList.currentText() == "Gaussian":
                        params = np.array([self.peak1_paramsEdits[0].value(), self.peak1_paramsEdits[1].value(), self.peak1_paramsEdits[2].value(), self.peak1_paramsEdits[3].value(), self.peak2_paramsEdits[0].value(), self.peak2_paramsEdits[1].value(), self.peak2_paramsEdits[2].value() ])
                        self.plotWidget.plot(self.x_fit, lineshape.Double_Gaussian(params, self.x_fit), pen=self.estPen)
                    else:
                        params = np.array([self.peak1_paramsEdits[0].value(), self.peak1_paramsEdits[1].value(), self.peak1_paramsEdits[2].value(), self.peak1_paramsEdits[3].value(), self.peak2_paramsEdits[0].value(), self.peak2_paramsEdits[1].value(), self.peak2_paramsEdits[2].value(), 0.9, 4 ])
                        self.plotWidget.plot(self.x_fit, lineshape.Double_Crystalball(params, self.x_fit), pen=self.estPen)


    def fitSingleGaussianPeak(self):
        params = np.array([self.peak1_paramsEdits[0].value(), self.peak1_paramsEdits[1].value(), self.peak1_paramsEdits[2].value(), self.peak1_paramsEdits[3].value() ])

        model = odr.Model(lineshape.Single_Gaussian)
        odr_data = odr.RealData(self.x_cut, self.y_cut )#, sy = np.sqrt(self.y_cut))
        myodr = odr.ODR(odr_data, model, beta0 = params)
        myodr.set_job(fit_type = 0)
        myoutput = myodr.run()

        self.params = myoutput.beta
        params_cov  = myoutput.cov_beta

        params_output = []
        params_red_chi2 = myoutput.res_var

        content, error = integrate.quad(lineshape.Single_Gaussian_integrand, self.x_fit[0], self.x_fit[-1], args=tuple(self.params))

        params_output.append(content)
        params_output.append(params_red_chi2)
        self.params_outputs = np.array(params_output)

        for i in range(len(self.peak1_paramsEdits)):
            self.peak1_paramsEdits[i].setValue(self.params[i])
        self.peak1_contentEdit.setText(str(self.params_outputs[0]))
        self.peak1_chi2Edit.setText(str(self.params_outputs[1]))

        self.plotWidget.plot(clear=True)
        self.plotWidget.plot(self.x,self.y[:-1],stepMode=True)
        fitPen = pg.mkPen(color=(204, 0, 0), width = 2)
        self.plotWidget.addItem(self.vLine, ignoreBounds=True)
        self.plotWidget.addItem(self.hLine, ignoreBounds=True)
        self.plotWidget.plot(self.x_fit, lineshape.Single_Gaussian(self.params, self.x_fit), pen=fitPen)



    def fitDoubleGaussianPeak(self):
        params = np.array([self.peak1_paramsEdits[0].value(), self.peak1_paramsEdits[1].value(), self.peak1_paramsEdits[2].value(), self.peak1_paramsEdits[3].value(), self.peak2_paramsEdits[0].value(), self.peak2_paramsEdits[1].value(), self.peak2_paramsEdits[2].value() ])

        model = odr.Model(lineshape.Double_Gaussian)
        odr_data = odr.RealData(self.x_cut, self.y_cut)#, sy = np.sqrt(self.y_cut))
        myodr = odr.ODR(odr_data, model, beta0 = params)
        myodr.set_job(fit_type = 0)
        myoutput = myodr.run()

        self.params = myoutput.beta
        params_cov  = myoutput.cov_beta

        params_output = []
        params_red_chi2 = myoutput.res_var

        params_single = self.params[:4]
        index = [4,5,6,3]
        params_double = self.params[index]

        content1, error1 = integrate.quad(lineshape.Single_Gaussian_integrand, self.x_fit[0], self.x_fit[-1], args=tuple(params_single))
        content2, error2 = integrate.quad(lineshape.Single_Gaussian_integrand, self.x_fit[0], self.x_fit[-1], args=tuple(params_double))

        params_output.append(content1)
        params_output.append(params_red_chi2)
        params_output.append(content2)
        self.params_outputs = np.array(params_output)

        for i in range(len(self.params)):
            if i < 4:
                self.peak1_paramsEdits[i].setValue(self.params[i])
            if 3 < i < 7:
                self.peak2_paramsEdits[i-4].setValue(self.params[i])

        self.peak1_contentEdit.setText(str(self.params_outputs[0]))
        self.peak1_chi2Edit.setText(str(self.params_outputs[1]))
        self.peak2_contentEdit.setText(str(self.params_outputs[2]))

        self.plotWidget.plot(clear=True)
        histo = pg.PlotCurveItem(self.x,self.y[:-1],stepMode=True)
        self.plotWidget.addItem(histo)

        if self.PeaksButton.isChecked() == False:
            fitPen = pg.mkPen(color=(204, 0, 0), width = 2)
            self.plotWidget.plot(self.x_fit, lineshape.Double_Gaussian(self.params, self.x_fit), pen=fitPen)
        if self.PeaksButton.isChecked() == True:
            fitPen = pg.mkPen(color=(204, 0, 0), width = 2, style=QtCore.Qt.DashLine)
            params1 = np.array([self.peak1_paramsEdits[0].value(), self.peak1_paramsEdits[1].value(), self.peak1_paramsEdits[2].value(), self.peak1_paramsEdits[3].value() ])
            params2 = np.array([self.peak2_paramsEdits[0].value(), self.peak2_paramsEdits[1].value(), self.peak2_paramsEdits[2].value(), self.peak1_paramsEdits[3].value() ])
            self.plotWidget.plot(self.x_fit, lineshape.Single_Gaussian(params1, self.x_fit), pen=fitPen)
            self.plotWidget.plot(self.x_fit, lineshape.Single_Gaussian(params2, self.x_fit), pen=fitPen)
        self.plotWidget.addItem(self.vLine, ignoreBounds=True)
        self.plotWidget.addItem(self.hLine, ignoreBounds=True)


    def fitSingleCrystalBallPeak(self):
        params = np.array([self.peak1_paramsEdits[0].value(), self.peak1_paramsEdits[1].value(), self.peak1_paramsEdits[2].value(), self.peak1_paramsEdits[3].value(), 0.9, 4 ])

        model = odr.Model(lineshape.Single_Crystalball)
        odr_data = odr.RealData(self.x_cut, self.y_cut)#, sy = np.sqrt(self.y_cut))
        myodr = odr.ODR(odr_data, model, beta0 = params)
        myodr.set_job(fit_type = 0)
        myoutput = myodr.run()

        self.params = myoutput.beta
        params_cov  = myoutput.cov_beta

        params_output = []
        params_red_chi2 = myoutput.res_var

        content, error = integrate.quad(lineshape.Single_Crystalball_integrand, self.x_fit[0], self.x_fit[-1], args=tuple(self.params))

        params_output.append(content)
        params_output.append(params_red_chi2)
        self.params_outputs = np.array(params_output)

        for i in range(len(self.peak1_paramsEdits)):
            self.peak1_paramsEdits[i].setValue(self.params[i])
        self.peak1_contentEdit.setText(str(self.params_outputs[0]))
        self.peak1_chi2Edit.setText(str(self.params_outputs[1]))

        self.plotWidget.plot(clear=True)
        self.plotWidget.plot(self.x,self.y[:-1],stepMode=True)
        fitPen = pg.mkPen(color=(204, 0, 0), width = 2)
        self.plotWidget.plot(self.x_fit, lineshape.Single_Crystalball(self.params, self.x_fit), pen=fitPen)
        self.plotWidget.addItem(self.vLine, ignoreBounds=True)
        self.plotWidget.addItem(self.hLine, ignoreBounds=True)


    def fitDoubleCrystalBallPeak(self):
        params = np.array([self.peak1_paramsEdits[0].value(), self.peak1_paramsEdits[1].value(), self.peak1_paramsEdits[2].value(), self.peak1_paramsEdits[3].value(), self.peak2_paramsEdits[0].value(), self.peak2_paramsEdits[1].value(), self.peak2_paramsEdits[2].value(), 0.9, 4 ])

        model = odr.Model(lineshape.Double_Crystalball)
        odr_data = odr.RealData(self.x_cut, self.y_cut)#, sy = np.sqrt(self.y_cut))
        myodr = odr.ODR(odr_data, model, beta0 = params)
        myodr.set_job(fit_type = 0)
        myoutput = myodr.run()

        self.params = myoutput.beta
        params_cov  = myoutput.cov_beta

        params_output = []
        params_red_chi2 = myoutput.res_var

        single_index = [0,1,2,3,7,8]
        params_single = self.params[single_index]
        double_index = [4,5,6,3,7,8]
        params_double = self.params[double_index]

        content1, error1 = integrate.quad(lineshape.Single_Crystalball_integrand, self.x_fit[0], self.x_fit[-1], args=tuple(params_single))
        content2, error2 = integrate.quad(lineshape.Single_Crystalball_integrand, self.x_fit[0], self.x_fit[-1], args=tuple(params_double))

        params_output.append(content1)
        params_output.append(params_red_chi2)
        params_output.append(content2)
        self.params_outputs = np.array(params_output)

        for i in range(len(self.params)):
            if i < 4:
                self.peak1_paramsEdits[i].setValue(self.params[i])
            if 3 < i < 7:
                self.peak2_paramsEdits[i-4].setValue(self.params[i])

        self.peak1_contentEdit.setText(str(self.params_outputs[0]))
        self.peak1_chi2Edit.setText(str(self.params_outputs[1]))
        self.peak2_contentEdit.setText(str(self.params_outputs[2]))

        self.plotWidget.plot(clear=True)
        histo = pg.PlotCurveItem(self.x,self.y[:-1],stepMode=True)
        self.plotWidget.addItem(histo)

        if self.PeaksButton.isChecked() == False:
            fitPen = pg.mkPen(color=(204, 0, 0), width = 2)
            self.plotWidget.plot(self.x_fit, lineshape.Double_Crystalball(self.params, self.x_fit), pen=fitPen)
        if self.PeaksButton.isChecked() == True:
            fitPen = pg.mkPen(color=(204, 0, 0), width = 2, style=QtCore.Qt.DashLine)
            params1 = np.array([self.peak1_paramsEdits[0].value(), self.peak1_paramsEdits[1].value(), self.peak1_paramsEdits[2].value(), self.peak1_paramsEdits[3].value(), 0.9, 4 ])
            params2 = np.array([self.peak2_paramsEdits[0].value(), self.peak2_paramsEdits[1].value(), self.peak2_paramsEdits[2].value(), self.peak1_paramsEdits[3].value(), 0.9, 4 ])
            self.plotWidget.plot(self.x_fit, lineshape.Single_Crystalball(params1, self.x_fit), pen=fitPen)
            self.plotWidget.plot(self.x_fit, lineshape.Single_Crystalball(params2, self.x_fit), pen=fitPen)
        self.plotWidget.addItem(self.vLine, ignoreBounds=True)
        self.plotWidget.addItem(self.hLine, ignoreBounds=True)

    def calcRate(self):
        if len(self.timestamps) != 0:
            TotalTime = np.sort(self.timestamps)[-1]
            self.TotalTimeEdit.setText(str(TotalTime))

            TotalCounts = 0
            for i in range(len(self.y)):
                TotalCounts += self.y[i]
            self.TotalCountsEdit.setText(str(TotalCounts))

            TotalRate = float(TotalCounts)/float(TotalTime)
            self.CountRateEdit.setText(str(TotalRate))
        else:
            print("No events for detector "+self.DetTypeList.currentText()+" in data file")

    def calcGatedRate(self):

        self.plotWidget.addItem(self.vLine, ignoreBounds=True)
        self.plotWidget.addItem(self.hLine, ignoreBounds=True)

        mask = np.bitwise_and(float(self.Gate_iEdit.value()) <= self.x, self.x <= float(self.Gate_fEdit.value()))
        y_gated = self.y[mask]

        GatedCounts = 0
        for i in range(len(y_gated)):
            GatedCounts += y_gated[i]
        self.GatedCountsEdit.setText(str(GatedCounts))

        TotalTime = np.sort(self.timestamps)[-1]

        GatedRate = float(GatedCounts)/float(TotalTime)
        self.GatedRateEdit.setText(str(GatedRate))

    def findCoincs(self):
                                            # Energy        Timestamp       Detector
        self.coinc_data = np.column_stack((self.data[:,3], self.data[:,1], self.data[:,2]))
        self.coinc_data = np.array(self.coinc_data, dtype=float)
        data = self.coinc_data[self.coinc_data[:,1].argsort()] # Tiem sort the data
        self.time_delta = data[1:,1] - data[:-1,1] # Calcualte delta-t between sequential events
        # Create coincidence array [Energy, Detector, Energy, Detector, Time_delta]
        self.coincidences = np.column_stack((data[:-1,0], data[:-1,2], data[1:,0], data[1:,2], self.time_delta))

    def recalcCoincs(self):


        time_start  = self.TimeStartEdit.value()
        time_stop   = self.TimeStopEdit.value()
        self.det_plot = float(self.Dets[self.DetectorPlotList.currentIndex()])
        self.det_gate = float(self.Dets[self.DetectorGateList.currentIndex()])
        self.e_gate_i = self.EnergyGate_iEdit.value()
        self.e_gate_f = self.EnergyGate_fEdit.value()

        self.coinc      = self.coincidences[(self.coincidences[:,4] >= time_start) & (self.coincidences[:,4] <= time_stop) & (self.coincidences[:,1] == self.det_plot) & (self.coincidences[:,3] == self.det_gate) & (self.coincidences[:,2] >= self.e_gate_i) & (self.coincidences[:,2] <= self.e_gate_f)]

        y_coinc, x  = np.histogram(self.coinc[:,0], self.NoBinsPlotEdit.value())
        x_coinc = x[0] + np.diff(x).cumsum()

        self.plotWidget.plot(clear=True)
        self.plotWidget.plot(x_coinc, y_coinc[:-1], stepMode=True, pen='r')
        self.plotWidget.addItem(self.vLine, ignoreBounds=True)
        self.plotWidget.addItem(self.hLine, ignoreBounds=True)

        self.EventsCoincEdit.setText(str(len(self.coinc)))


    def show2Dhist(self):

        H, xedges, yedges = np.histogram2d(self.coinc[:,0], self.coinc[:,2], [self.NoBinsPlotEdit.value(), self.NoBinsGateEdit.value()])
        H = np.rot90(H)
        H = np.flipud(H)
        Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero

        plt.figure()
        plt.pcolormesh(xedges,yedges,Hmasked)
        plt.xlabel(self.DetectorPlotList.currentText()+" energy (keV)")
        plt.ylabel(self.DetectorGateList.currentText()+" energy (keV)")
        plt.colorbar()
        plt.title(self.DetectorPlotList.currentText()+" - "+self.DetectorGateList.currentText()+" coincidences")
        plt.show()

    def showdthist(self):

        coinc_dt = self.coincidences[(self.coincidences[:,1] == self.det_plot) & (self.coincidences[:,3] == self.det_gate) & (self.coincidences[:,2] >= self.e_gate_i) & (self.coincidences[:,2] <= self.e_gate_f)]

        plt.hist(coinc_dt[:,4], bins=5000, histtype='step', color='b')
        plt.xlabel(r'Time between successive events ($\Delta t$)')
        plt.ylabel("Counts per bin")
        plt.show()


    def mouseMoved(self, evt):
        mousePoint = self.plotWidget.plotItem.vb.mapSceneToView(evt[0])
        self.label.setText("x = %.1f, y = %.1f" % (mousePoint.x(), mousePoint.y()))
        self.vLine.setPos(mousePoint.x()), self.hLine.setPos(mousePoint.y())

def main():

    app = QtGui.QApplication(sys.argv)
    ex = FitDSSData()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
