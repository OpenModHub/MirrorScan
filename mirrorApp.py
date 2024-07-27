import sys
import yaml
from PySide6.QtWidgets import QApplication, QFileDialog, QLabel, QVBoxLayout, QWidget, QProgressBar, QMessageBox
from PySide6.QtCore import QObject, QThread, Signal, Slot
from PySide6.QtGui import QTransform
import pyqtgraph.opengl as gl
import pyqtgraph as pg
from pyqtgraph import functions as fn
import numpy as np
import os
import asyncio
from time import sleep
import datetime
from timeit import default_timer as timer

            
try:
    import nea_tools
    offline_mode = False
except:
    # TODO replace with logger
    print("nea_tools module not found, working offline")
    offline_mode = True

current_folder = os.getcwd()
ui_file = os.path.join(current_folder,'mirrorApp.ui')

uiclass, baseclass = pg.Qt.loadUiType(ui_file)

######## QT WORKING THREAD CLASS ############
class Worker(QObject):
    progress = Signal(int)
    advanced_progress = Signal(int)
    Zcompleted = Signal()
    completed = Signal()
    advanced_completed = Signal()
    started = Signal()
    advanced_started = Signal()
    status_update = Signal(str)

    def __init__(self):
        super().__init__()
        self.scan_map = mirror_scan()
        self.advanced_scan_map = advanced_scan()
        self.advanced_positions = []
        self.nea = None
        self.context = None
        self.motors = None
        self.Vector3D = None
        self.advanced_aborted = False

    @Slot()
    def do_scan(self):
        self.started.emit()
        # Create motor object
        p = self.motors.Mirror()
        if not p.is_active:
            p.activate()
        # Set sampling interval
        self.context.Microscope.Py.SetSamplingTime(50)
        # Set motor speed
        safe_v = self.context.Microscope.Py.MirrorMotorVelocityInContacting
        v = self.Vector3D(safe_v,safe_v,safe_v)
        self.context.Microscope.Py.SetActiveMotorVelocityXyz(v)
        # Update current position
        self.context.Microscope.RefreshActiveMotorPositionXyzAsync().Wait()
        current_pos = p.absolute_position
        self.scan_map.center_point = current_pos
        self.status_update.emit(f'Mirror position BEFORE movement: {current_pos}')

        # Calculate mirror coordinates for movement
        xs = current_pos[0] + np.linspace(-self.scan_map.sizeX/2,self.scan_map.sizeX/2,self.scan_map.Nx)
        ys = current_pos[1] + np.linspace(-self.scan_map.sizeY/2,self.scan_map.sizeY/2,self.scan_map.Ny)

        if self.scan_map.Nz == 1:
            # self.scan_map.Nz = 1
            zs = np.array([current_pos[2]])
        else:
            zs = current_pos[2] + np.linspace(-self.scan_map.sizeZ/2,self.scan_map.sizeZ/2,self.scan_map.Nz)

        # SCANNING LOOP
        counter = 0
        counterZ = 0
        startime = timer()
        for idz, z in enumerate(zs):
            for idy, y in enumerate(ys):
                for idx, x in enumerate(xs):
                    counter += 1
                    self.context.Microscope.RefreshActiveMotorPositionXyzAsync().Wait()
                    posx = p.absolute_position[0]
                    posy = p.absolute_position[1]
                    posz = p.absolute_position[2]
                    dx = x-posx
                    dy = y-posy
                    dz = z-posz
                    p.go_relative(dx,dy,dz)
                    p.await_movement()
                    # Read optical channels
                    self.scan_map.O1A[idz,idy,idx] = self.context.Microscope.Py.OpticalAmplitude[1]
                    self.scan_map.O2A[idz,idy,idx] = self.context.Microscope.Py.OpticalAmplitude[2]
                    self.scan_map.O3A[idz,idy,idx] = self.context.Microscope.Py.OpticalAmplitude[3]
                    self.scan_map.O4A[idz,idy,idx] = self.context.Microscope.Py.OpticalAmplitude[4]
                    # Update real position
                    self.context.Microscope.RefreshActiveMotorPositionXyzAsync().Wait()
                    newx = p.absolute_position[0]
                    newy = p.absolute_position[1]
                    newz = p.absolute_position[2]
                    self.scan_map.X[idz,idy,idx] = newx
                    self.scan_map.Y[idz,idy,idx] = newy
                    self.scan_map.Z[idz,idy,idx] = newz
                    steptime = timer()
                    remtime = (steptime-startime)/counter*(self.scan_map.Nx*self.scan_map.Ny*self.scan_map.Nz-counter)
                    self.progress.emit(counterZ)
                    self.status_update.emit(f'X: {newx}, Y: {newy}, Z: {newz} Remaining time: {datetime.timedelta(seconds=remtime)}')
            counterZ += 1
            if self.scan_map.Nz > 1:
                self.Zcompleted.emit()
        sleep(0.5)

        # Go back to the original position
        self.context.Microscope.RefreshActiveMotorPositionXyzAsync().Wait()
        dx = current_pos[0]-p.absolute_position[0]
        dy = current_pos[1]-p.absolute_position[1]
        dz = current_pos[2]-p.absolute_position[2]
        p.go_relative(dx,dy,dz)
        p.await_movement()
        # Check position after going back:
        self.context.Microscope.RefreshActiveMotorPositionXyzAsync().Wait()
        current_pos = p.absolute_position
        self.scan_map.center_point = current_pos
        self.status_update.emit(f'Mirror position AFTER movement: {current_pos}')
        self.completed.emit()

    def go_to_center(self,map,p):
        # REAL FUNCTION
        # self.context.Microscope.RefreshActiveMotorPositionXyzAsync().Wait()
        # dx = map.center_point[0]-p.absolute_position[0]
        # dy = map.center_point[1]-p.absolute_position[1]
        # dz = map.center_point[2]-p.absolute_position[2]
        # p.go_relative(dx,dy,dz)
        # p.await_movement()
        # self.context.Microscope.RefreshActiveMotorPositionXyzAsync().Wait()
        # current_pos = p.absolute_position
        # self.map.center_point = current_pos
        # self.status_update.emit(f'Mirror position after RECENTERING: {self.map.center_point}')

        # FOR TEST
        map.center_point = [0.0,0.0,0.0]
        self.status_update.emit(f'Mirror position after RECENTERING: {map.center_point}')
        sleep(0.5)

    @Slot()
    def do_scan_advanced(self):
        self.advanced_started.emit()
        # Create motor object
        # p = self.motors.Mirror()
        # if not p.is_active:
        #     p.activate()
        # Set sampling interval
        # self.context.Microscope.Py.SetSamplingTime(50)
        # Set motor speed
        # safe_v = self.context.Microscope.Py.MirrorMotorVelocityInContacting
        # v = self.Vector3D(safe_v,safe_v,safe_v)
        # self.context.Microscope.Py.SetActiveMotorVelocityXyz(v)
        # Update current position
        # self.context.Microscope.RefreshActiveMotorPositionXyzAsync().Wait()
        # current_pos = p.absolute_position
        xpositions = self.advanced_positions[:,0]
        ypositions = self.advanced_positions[:,1]
        zpositions = self.advanced_positions[:,2]
        
        current_pos = [0.0,0.0,0.0]
        self.advanced_scan_map.center_point = current_pos
        self.status_update.emit(f'Mirror position BEFORE movement: {current_pos}')

        # Calculate mirror coordinates for movement
        xs = current_pos[0] + xpositions
        ys = current_pos[1] + ypositions
        zs = current_pos[2] + zpositions

        # SCANNING LOOP
        counter = 0
        startime = timer()
        for i in range(self.advanced_scan_map.Npoints):
            counter += 1
            # self.context.Microscope.RefreshActiveMotorPositionXyzAsync().Wait()
            # posx = p.absolute_position[0]
            # posy = p.absolute_position[1]
            # posz = p.absolute_position[2]
            if i == 0:
                dx = xs[i] - current_pos[0]
                dy = xs[i] - current_pos[1]
                dz = ys[i] - current_pos[2]
            else:
                dx = xs[i] - xs[i-1]
                dy = ys[i] - ys[i-1]
                dz = zs[i] - zs[i-1]
            
            # dx = xs[i]-posx
            # dy = ys[i]-posy
            # dz = zs[i]-posz
            # p.go_relative(dx,dy,dz)
            # p.await_movement()
            # print(f"Move relative: {dx},{dy},{dz}")

            # Read optical channels
            # self.advanced_scan_map.O1A[i] = self.context.Microscope.Py.OpticalAmplitude[1]
            # self.advanced_scan_map.O2A[i] = self.context.Microscope.Py.OpticalAmplitude[2]
            # self.advanced_scan_map.O3A[i] = self.context.Microscope.Py.OpticalAmplitude[3]
            # self.advanced_scan_map.O4A[i] = self.context.Microscope.Py.OpticalAmplitude[4]
            vv = np.array([xs[i],ys[i],zs[i]])
            self.advanced_scan_map.O1A.append(np.linalg.norm(vv))
            self.advanced_scan_map.O2A.append(np.linalg.norm(vv))
            self.advanced_scan_map.O3A.append(np.linalg.norm(vv))
            self.advanced_scan_map.O4A.append(np.linalg.norm(vv))

            # Update real position
            # self.context.Microscope.RefreshActiveMotorPositionXyzAsync().Wait()
            # newx = p.absolute_position[0]
            # newy = p.absolute_position[1]
            # newz = p.absolute_position[2]
            newx = xs[i]
            newy = ys[i]
            newz = zs[i]
            self.advanced_scan_map.X.append(newx)
            self.advanced_scan_map.Y.append(newy)
            self.advanced_scan_map.Z.append(newz)

            sleep(0.1)
            self.advanced_progress.emit(counter)
            steptime = timer()
            remtime = (steptime-startime)/counter*(self.advanced_scan_map.Npoints-counter)
            self.status_update.emit(f'X: {newx}, Y: {newy}, Z: {newz} Remaining time: {datetime.timedelta(seconds=remtime)}')
            if self.advanced_aborted:
                break

        sleep(0.5)
        # Go back to the original position
        # self.context.Microscope.RefreshActiveMotorPositionXyzAsync().Wait()
        # dx = current_pos[0]-p.absolute_position[0]
        # dy = current_pos[1]-p.absolute_position[1]
        # dz = current_pos[2]-p.absolute_position[2]
        # p.go_relative(dx,dy,dz)
        # p.await_movement()
        # Check position after going back:
        # self.context.Microscope.RefreshActiveMotorPositionXyzAsync().Wait()
        # current_pos = p.absolute_position
        # self.advanced_scan_map.center_point = current_pos
        if self.advanced_aborted:
            self.go_to_center(self.advanced_scan_map,p=None)
            self.status_update.emit(f'Scan was ABORTED!!! Tip position: {self.advanced_scan_map.center_point}')
            self.advanced_aborted = False
        else:
            self.status_update.emit(f'Mirror position AFTER movement: {current_pos}')
            self.advanced_completed.emit()

######## MAIN APPLICATION WINDOW CLASS ############
class MainWindow(uiclass, baseclass):
    work_requested = Signal()
    advanced_work_requested = Signal()

    pg.setConfigOptions(imageAxisOrder='row-major')

    def __init__(self):
        super().__init__()

        # Load UI
        self.setupUi(self)

        # Other attributes and flags
        self.offline_mode = offline_mode
        self.connected = False
        self.click_move_enabled = False
        # config parameters for SNOM connection
        self.config = None
        # To store the measured maps
        self.mirror_map = None
        self.advanced_map = None
        self.advanced_positions = None
        # To store the loaded maps
        self.loaded_map = None
        self.loaded_advanced_map = None
        # The map the is currently plotted
        self.advanced_map_to_plot = None
        # Flags
        self.sizes_linked = False
        self.step_sizes_linked = False
        self.currentZindex = 0
        
        # Stylize
        self.setWindowTitle('Focus scanner application')
        # self.setStyleSheet("background-color: white;")
        self.channel_comboBox.addItems(['O1A', 'O2A', 'O3A','O4A'])
        self.channel_comboBox.setCurrentText('O3A')
        self.advanced_channel_comboBox.addItems(['O1A', 'O2A', 'O3A','O4A'])
        self.advanced_channel_comboBox.setCurrentText('O3A')

        # Linking button label correction
        txt = "\U0001F517"
        self.linkSizeButton.setText(txt)
        self.linkStepSizeButton.setText(txt)
        self.linkButton_default_style_sheet = self.linkSizeButton.styleSheet()

        # Create the worker thread
        self.worker = Worker()
        self.worker_thread = QThread()
        # Connected worker thread signals
        self.worker.progress.connect(self.update_scan_progress)
        self.worker.completed.connect(self.scan_complete)
        self.worker.Zcompleted.connect(self.save_data(fname='temp.dat'))
        self.worker.advanced_progress.connect(self.update_advanced_scan_progress)
        self.worker.advanced_completed.connect(self.advanced_scan_complete)
        self.worker.advanced_started.connect(self.advanced_scan_started)
        self.worker.status_update.connect(self.status_bar_update)
        self.work_requested.connect(self.worker.do_scan)
        self.advanced_work_requested.connect(self.worker.do_scan_advanced)
        self.worker.moveToThread(self.worker_thread)
        self.worker_thread.start()
        
        # Create test data
        testdata = np.fromfunction(lambda i, j: (1+0.3*np.sin(i)) * (i)**2 + (j)**2, (100, 100))
        testdata = testdata * (1 + 0.2 * np.random.random(testdata.shape) )
        testdata = testdata.transpose()
        self.data_to_plot = testdata
        self.center_pos_rel = [0, 0]
        self.center_pos_abs = [50, 50]
        self.center_marker = [{'pos': [50, 50], 'data': 1}]

        # Create plot widget
        self.imItem = pg.ImageItem(image=testdata)                                                  # create an ImageItem
        self.plot_area.addItem(self.imItem)                                                         # add it to the PlotWidget
        self.cbar = self.plot_area.addColorBar(self.imItem, colorMap='CET-L9',rounding=0.01)        # Create a colorBarItem and add to the PlotWidget
        self.scatterItem = pg.ScatterPlotItem(size=15,
                                              pen=pg.mkPen(color=(255, 255, 255, 220), width=1.5),
                                              brush=pg.mkBrush(255, 255, 255, 120))
        self.scatterItem.addPoints(self.center_marker)
        self.plot_area.addItem(self.scatterItem)
        # self.plot_area.setBackground('w')
        self.plot_area.getAxis('left').setTextPen('black')
        self.plot_area.getAxis('bottom').setTextPen('black')
        self.plot_area.setBackground((200,200,200, 1))
        self.plot_area.setMouseTracking(True)                                                       # For cursor tracking
        self.imItem.hoverEvent = self.imageHoverEvent                                               # Attach event
        self.imItem.mouseClickEvent = self.imageClickEvent                                          # Attach event

        # Create the plot widget for the advanced mode with pyqtgraph OpenGL scatterplot
        self.gridItem = gl.GLGridItem()
        self.plot_area_advanced.addItem(self.gridItem)
        # Create some test data
        pos = np.random.uniform(low=-25, high=25, size=(100000,3))
        pos[0] = (0,0,0)

        d = []
        for i in range(np.shape(pos)[0]):
            d.append(1-np.linalg.norm(pos[i,:]))
        d = np.array(d)

        colors = np.ones((pos.shape[0], 4))
        minval = np.min(d)
        maxval = np.max(d)
        colors[:,0:3] = self.calculateColors(data = d, vmin = minval, vmax = maxval)
        alphas = (d-minval)/(maxval-minval)
        colors[:,3] = alphas

        self.plot3DItem = gl.GLScatterPlotItem(pos=pos, color=(1,1,1,1), size=0.1, pxMode=False)
        self.plot3DItem.setData(pos=pos, color=colors)
        self.plot_area_advanced.addItem(self.plot3DItem)
        # self.plot3DItem.setGLOptions('opaque')

        # Connect button signals
        self.choose_file_button.clicked.connect(self.choose_file)
        self.datascroll_spinBox.valueChanged.connect(self.data_scroll)
        self.channel_comboBox.currentIndexChanged.connect(self.channel_change)
        self.advanced_channel_comboBox.currentIndexChanged.connect(self.advanced_channel_change)
        self.linkSizeButton.clicked.connect(self.link_scan_size)
        self.linkStepSizeButton.clicked.connect(self.link_scan_step_size)
        self.load_coords_button.clicked.connect(self.load_advanced_points)
        self.load_meas_adv_button.clicked.connect(self.load_advanced_data)
        self.abort_adv_button.clicked.connect(self.abort_advanced)
        self.show_coords_button.clicked.connect(self.show_advanced_points)

        # Check if config file is modified
        self.check_config_file()

        if not offline_mode:
            self.scan_button.clicked.connect(self.start_scan)
            self.connect_snom_button.clicked.connect(self.connect_to_neasnom)
            self.move_to_button.clicked.connect(self.enable_move_to_point)
            self.save_button.clicked.connect(self.save_data)

        self.star_scan_adv_button.clicked.connect(self.start_advanced_scan)

        if offline_mode:
            # simaple buttons
            self.statusbar.showMessage(u"\u26A0 nea_tools module not found, running in display-only mode.")
            self.connect_snom_button.setEnabled(False)
            self.move_to_button.setEnabled(False)
            self.scan_button.setEnabled(False)
            self.save_button.setEnabled(False)
            # advanced buttons
            # self.star_scan_adv_button.setEnable(False)
            self.abort_adv_button.setEnabled(False)

    def check_config_file(self):   # load config
        with open('config.yaml', 'r') as file:
            self.config = yaml.safe_load(file)
            if (self.config['fingerprint'] == 'CHANGEMEE') or (self.config['path_to_dll'] == r"CHANGEMEE"):
                msg = QMessageBox()
                msg.setWindowTitle("Configuration missing")
                msg.setText("You have to set up neaSNOM configuration before use")
                msg.setIcon(QMessageBox.Critical)
                msg.setStandardButtons(QMessageBox.Ok|QMessageBox.Cancel)
                buttonConnect = msg.button(QMessageBox.Ok)
                buttonConnect.setText('Ok')
                msg.setInformativeText("Click 'Ok' and set the parameters in the config.yaml file or click 'Cancel' to continue in offline mode")
                button = msg.exec()
                if button == QMessageBox.Ok:
                    sys.exec()
                elif button == QMessageBox.Cancel:
                    self.offline_mode = True

    def connect_to_neasnom(self):
        if "nea_tools" not in sys.modules:
            return

        self.path_to_dll = ''# yaml.load('config.yaml')
        path_to_dll = self.config['path_to_dll']
        fingerprint = self.config['fingerprint']
        host = 'nea-server'
        if self.connected:
            print('\nDisconnecting from neaServer!')
            nea_tools.disconnect()
            self.connected = False
            self.connect_snom_button.setText("Connect to neaSNOM")
            self.statusbar.showMessage("Disconnected from SNOM")
        else:
            loop = asyncio.get_event_loop()
            try:
                loop.run_until_complete(nea_tools.connect(host, fingerprint, path_to_dll))
            except ConnectionError:
                print("Could not connect!")
            try:
                from neaspec import context
                import Nea.Client.SharedDefinitions as nea
                from nea_tools.microscope import motors
                self.Vector3D = nea.Geometry.Vector3D
                self.Point3D = nea.Geometry.Point3D
            except ModuleNotFoundError:
                raise ConnectionError('Connection refused or timeout. Retry to connect again.')
            else:
                self.connected = True
                self.statusbar.showMessage("Connected to neaSNOM")
                self.connect_snom_button.setText("Disconnect neaSNOM")

            self.context = context
            self.nea = nea
            self.motors = motors
            return context, nea

    def choose_file(self):
        fname = QFileDialog.getOpenFileName(self, "Choose file","","Datatext files (*.txt *.dat)")
        self.file_name = fname[0]
        try:
            self.load_data()
            self.statusbar.showMessage(f'{fname} is loaded')
        except:
            self.statusbar.showMessage(f'No file was loaded')
    
    def load_data(self):
        # Create scan object
        self.loaded_map = mirror_scan()
        # Read header lines
        nlines = 6
        with open(self.file_name, 'r') as file: header_lines = [file.readline().strip() for _ in range(nlines)]

        for header_line in header_lines:
            idx = header_line.find("=")
            text = header_line[2:idx-1]
            number = float(header_line[idx+2:])
            if text == 'SizeX':
                self.loaded_map.sizeX = number
            elif text == 'SizeY':
                self.loaded_map.sizeY = number
            elif text == 'SizeZ':
                self.loaded_map.sizeZ = number
            elif text == 'StepX':
                self.loaded_map.step_sizeX = number
            elif text == 'StepY':
                self.loaded_map.step_sizeY = number
            elif text == 'StepZ':
                self.loaded_map.step_sizeZ = number

        self.loaded_map.recalc_size()
        self.Zaxis = np.linspace(-self.loaded_map.sizeZ/2,self.loaded_map.sizeZ/2,self.loaded_map.Nz)

        # Load data section of the file and reshape it to the right size
        data = np.loadtxt(self.file_name)
        self.loaded_map.X = np.reshape(data[:,0],(self.loaded_map.Nz,self.loaded_map.Nx,self.loaded_map.Ny))
        self.loaded_map.Y = np.reshape(data[:,1],(self.loaded_map.Nz,self.loaded_map.Nx,self.loaded_map.Ny))
        self.loaded_map.Z = np.reshape(data[:,2],(self.loaded_map.Nz,self.loaded_map.Nx,self.loaded_map.Ny))
        self.loaded_map.O1A = np.reshape(data[:,3],(self.loaded_map.Nz,self.loaded_map.Nx,self.loaded_map.Ny))
        self.loaded_map.O2A = np.reshape(data[:,4],(self.loaded_map.Nz,self.loaded_map.Nx,self.loaded_map.Ny))
        self.loaded_map.O3A = np.reshape(data[:,5],(self.loaded_map.Nz,self.loaded_map.Nx,self.loaded_map.Ny))
        self.loaded_map.O4A = np.reshape(data[:,6],(self.loaded_map.Nz,self.loaded_map.Nx,self.loaded_map.Ny))
        # Create map object for the loaded data

        self.center_pos_rel = [0,0]
        self.center_marker = [{'pos': self.center_pos_rel, 'data': 1}]

        try:
            self.set_display_data(self.loaded_map)
            self.update_image()
        except:
            pass
        else:
            self.mirror_map = None

    def load_advanced_data(self):
        fname = QFileDialog.getOpenFileName(self, "Choose file","","Datatext files (*.txt *.dat)")
        file_name = fname[0]
        # Create scan object
        try:
            self.loaded_advanced_map = advanced_scan()
            # Load data
            data = np.loadtxt(file_name)
            self.loaded_advanced_map.X = data[:,0]
            self.loaded_advanced_map.Y = data[:,1]
            self.loaded_advanced_map.Z = data[:,2]
            self.loaded_advanced_map.O1A = data[:,3]
            self.loaded_advanced_map.O2A = data[:,4]
            self.loaded_advanced_map.O3A = data[:,5]
            self.loaded_advanced_map.O4A = data[:,6]

            self.set_advanced_display(self.loaded_advanced_map)
            self.update_advanced_plot()
        except IOError as e:
            self.loaded_advanced_map = None
            self.status_bar_update("No file was loaded")

    def load_advanced_points(self):
        fname = QFileDialog.getOpenFileName(self, "Choose file","","Datatext files (*.txt *.dat)")
        file_name = fname[0]
        try:
            self.advanced_positions = np.loadtxt(file_name)
            self.advanced_positions *= 1000
            print(f"Number of positions: {len(self.advanced_positions)}")
            self.number_of_points_label.setText(f"Number of positions: {len(self.advanced_positions)}")
            self.xrange_label.setText(f"X range: {np.min(self.advanced_positions[:,0])} - {np.max(self.advanced_positions[:,0])}")
            self.yrange_label.setText(f"Y range: {np.min(self.advanced_positions[:,1])} - {np.max(self.advanced_positions[:,1])}")
            self.zrange_label.setText(f"Z range: {np.min(self.advanced_positions[:,2])} - {np.max(self.advanced_positions[:,2])}")
        except IOError as e:
            self.status_bar_update("No file was loaded")

    def set_advanced_display(self,map):
        self.advanced_map_to_plot = map

    def update_advanced_plot(self):
        pos = np.zeros((len(self.advanced_map_to_plot.X), 3))
        # color = np.ones((len(self.advanced_map_to_plot.X), 4))
        pos[:,0] = np.asarray(self.advanced_map_to_plot.X)/1000
        pos[:,1] = np.asarray(self.advanced_map_to_plot.Y)/1000
        pos[:,2] = np.asarray(self.advanced_map_to_plot.Z)/1000
        colors = np.ones((pos.shape[0], 4))

        to_plot = 1 - np.array(getattr(self.advanced_map_to_plot,self.advanced_channel_comboBox.currentText()))
        # to_plot = 1-np.array(self.advanced_map_to_plot.O2A)

        minval = np.min(to_plot)
        maxval = np.max(to_plot)

        colors[:,0:3] = self.calculateColors(data = to_plot, vmin = minval, vmax = maxval)
        alphas = (to_plot-minval)/(maxval-minval)
        colors[:,3] = alphas
        sizescaling = 1

        self.plot3DItem.setData(pos=pos, color=colors, size=alphas*sizescaling)

    def calculateColors(self, data, vmin, vmax, colormapname = "viridis", n=256):
        # Make lookuptable
        cm = pg.colormap.get(colormapname)
        lut = cm.getLookupTable(nPts=n, alpha=False)
        # Scale values between 0-255
        normalized_values = (data-vmin)/(vmax-vmin)
        scaled_values = (normalized_values*255).astype(int)
        scaled_values = np.clip(scaled_values, 0, 255)
        
        rgba_colors = lut[scaled_values]
        rgba_colors = rgba_colors.astype(float)/255.0

        return rgba_colors
    
    def scaleTo255(self, data, vmin, vmax):
        normalized_values = (data-vmin)/(vmax-vmin)
        scaled_values = (normalized_values*255).astype(int)
        scaled_values = np.clip(scaled_values, 0, 255)

        return scaled_values

    def show_advanced_points(self):
        if self.advanced_positions is None:
            self.status_bar_update('No positions were loaded for advanced scan!')
            print("No positions loaded")
            msg = QMessageBox(self)
            msg.setWindowTitle("Missing point positions")
            msg.setText("Positions of measurement points are not defined")
            msg.setIcon(QMessageBox.Critical)
            msg.setStandardButtons(QMessageBox.Ok|QMessageBox.Cancel)
            buttonConnect = msg.button(QMessageBox.Ok)
            buttonConnect.setText('Choose file')
            msg.setInformativeText("Click 'Choose file' to browse for text file defining point positions")
            button = msg.exec()
            if button == QMessageBox.Ok:
                self.load_advanced_points()
            else:
                return None
        
        point_map = advanced_scan()
        # Load data
        point_map.X = self.advanced_positions[:,0]
        point_map.Y = self.advanced_positions[:,1]
        point_map.Z = self.advanced_positions[:,2]
        for i in range(np.shape(point_map.X)[0]):
            vv = np.array([point_map.X[i],point_map.Y[i],point_map.Z[i]])
            point_map.O1A.append(np.linalg.norm(vv))
            point_map.O2A.append(pow(np.linalg.norm(vv),2))
            point_map.O3A.append(np.sqrt(np.linalg.norm(vv)))
            point_map.O4A.append(np.log(np.linalg.norm(vv)))
        self.set_advanced_display(point_map)
        self.update_advanced_plot()

    def set_display_data(self,map):
        self.meas_data = np.array(getattr(map,self.channel_comboBox.currentText()))
        # self.meas_data = self.meas_data.reshape((map.Nz,map.Nx,map.Ny))
        if map.Nz == 1:
            self.data_to_plot = self.meas_data[0,:,:]
            self.datascroll_spinBox.setValue(0)
            self.datascroll_spinBox.setRange(0,0)
            self.datascroll_spinBox.setEnabled(False)
            self.cbar.setLevels(values = (np.min(self.data_to_plot),np.max(self.data_to_plot)))
        else:
            self.zeroZ_data = self.meas_data[self.currentZindex,:,:]
            self.data_to_plot = self.zeroZ_data
            self.datascroll_spinBox.setValue(self.currentZindex) #Here should come counterZ
            self.datascroll_spinBox.setRange(0, self.currentZindex)
            self.datascroll_spinBox.setEnabled(True)
            self.cbar.setLevels(values = (np.min(self.zeroZ_data),np.max(self.zeroZ_data)))

        #Set up the axis values by transforming the image
        tr = QTransform()                                                               # prepare ImageItem transformation:
        tr.translate(-map.sizeX/2/1000,-map.sizeY/2/1000)                                         # move 3x3 image to locate center at axis origin
        tr.scale(map.sizeX/map.Nx/1000, map.sizeY/map.Ny/1000)                                    # scale horizontal and vertical axes
        self.imItem.setTransform(tr)
        self.plot_area.getAxis('bottom').setLabel('X position / μm')
        self.plot_area.getAxis('left').setLabel('Y position / μm')
        self.plot_area.showAxes(True)
        self.plot_area.setAspectLocked(True)

    def update_image(self):
        lastlevels = self.cbar.levels()
        self.imItem.setImage(image = self.data_to_plot)
        self.cbar.setLevels(values = lastlevels)
        self.scatterItem.clear()
        self.scatterItem.addPoints(self.center_marker)
        self.click_move_enabled = False
        self.move_to_button.setEnabled(True)

    def data_scroll(self):
        if (self.loaded_map == None) & (self.mirror_map == None):
            print(f'Error')
        else:
            index = self.datascroll_spinBox.value()
            print(f'Z index: {index}, map size: {np.size(self.meas_data[index,:,:])}')
            self.data_to_plot = self.meas_data[index,:,:]
            self.update_image()
            # self.Zplane_label.setText(f"Displayed Z plane: {self.Zaxis[index]} nm")
    
    def channel_change(self):
        if (self.mirror_map == None) & (self.loaded_map == None):
            self.statusbar.showMessage(f"No scan data to display!")
        elif self.loaded_map == None:
            self.set_display_data(self.mirror_map)
            self.update_image()
            self.statusbar.showMessage(f"Channel changed to {self.channel_comboBox.currentText()}")
        else:
            self.set_display_data(self.loaded_map)
            self.update_image()
            self.statusbar.showMessage(f"Channel changed to {self.channel_comboBox.currentText()}")

    def advanced_channel_change(self):
        if (self.advanced_map_to_plot == None):
            self.statusbar.showMessage(f"No scan data to display!")
        else:
            self.update_advanced_plot()
            self.statusbar.showMessage(f"Channel changed to {self.channel_comboBox.currentText()}")

    def imageHoverEvent(self,event):
        # Show the position, pixel, and value under the mouse cursor.
        if event.isExit():
            self.plot_area.setTitle("")
            return
        pos = event.pos()
        i, j = pos.y(), pos.x()
        i = int(np.clip(i, 0, self.data_to_plot.shape[0] - 1))
        j = int(np.clip(j, 0, self.data_to_plot.shape[1] - 1))
        val = self.data_to_plot[i, j]
        ppos = self.imItem.mapToParent(pos)
        x, y = ppos.x(), ppos.y()
        self.plot_area.setTitle("pos: (%0.1f, %0.1f)  pixel: (%d, %d)  value: %.3g" % (x, y, i, j, val))

    def imageClickEvent(self, event):
        if self.click_move_enabled:
            # Create motor object
            p = self.motors.Mirror()
            if not p.is_active:
                p.activate()
            self.context.Microscope.RefreshActiveMotorPositionXyzAsync().Wait()
            pos_before = p.absolute_position
            self.status_bar_update(f"Absolute BEFORE center move: {pos_before}")
            # Get mouse position
            pos = event.pos()
            ppos = self.imItem.mapToParent(pos)
            x, y = ppos.x(), ppos.y()
            dx = x - self.center_pos_rel[0]
            dy = y - self.center_pos_rel[1]
            # Go to position
            p.go_relative(dx*1000,dy*1000,0)
            p.await_movement()
            self.status_bar_update(f"Relative move to {[x,y]}")
            # Check new position
            self.context.Microscope.RefreshActiveMotorPositionXyzAsync().Wait()
            self.center_pos_abs = p.absolute_position
            self.status_bar_update(f"Absolute AFTER center move: {self.center_pos_abs}")
            # Replace center marker
            realdx = (self.center_pos_abs[0] - pos_before[0])/1000
            realdy = (self.center_pos_abs[1] - pos_before[1])/1000
            self.center_pos_rel = [self.center_pos_rel[0] + realdx, self.center_pos_rel[1] + realdy]
            self.center_marker = [{'pos': self.center_pos_rel, 'data': 1}]
            self.scatterItem.clear()
            self.scatterItem.addPoints(self.center_marker)
            self.click_move_enabled = False
            self.move_to_button.setEnabled(True)
        else:
            pass

    def start_scan(self):
        # Create map object and set up scan parameters
        self.mirror_map = mirror_scan()
        self.mirror_map.step_sizeX = self.stepX_spinBox.value()*1000 #in nm
        self.mirror_map.step_sizeY = self.stepY_spinBox.value()*1000
        self.mirror_map.step_sizeZ = self.stepZ_spinBox.value()*1000
        self.mirror_map.sizeX = self.sizeX_spinBox.value()*1000
        self.mirror_map.sizeY = self.sizeY_spinBox.value()*1000
        self.mirror_map.sizeZ = self.sizeZ_spinBox.value()*1000
        self.mirror_map.recalc_size()
        self.mirror_map.create_array()
        self.Zaxis = np.linspace(-self.mirror_map.sizeZ/2,self.mirror_map.sizeZ/2,self.mirror_map.Nz)
        print(f'StepX: {self.mirror_map.step_sizeX} StepY: {self.mirror_map.step_sizeY}')
        # Send the map object to worker object
        self.worker.scan_map = self.mirror_map
        # Check if connected
        if self.connected:
            # Pass SDK objects to worker thread
            self.worker.nea = self.nea
            self.worker.context = self.context
            self.worker.motors = self.motors
            self.worker.Vector3D = self.Vector3D
            # Emit Signal to start scan at worker thread Slot
            self.work_requested.emit()
            self.connect_snom_button.setEnabled(False)
        else:
            self.status_bar_update('Connect to neaSNOM before scanning!')
            print("No device was found")
            msg = QMessageBox(self)
            msg.setWindowTitle("No connection!")
            msg.setText("Not connected to SNOM!")
            msg.setIcon(QMessageBox.Critical)
            msg.setStandardButtons(QMessageBox.Ok|QMessageBox.Cancel)
            buttonConnect = msg.button(QMessageBox.Ok)
            buttonConnect.setText('Connect')
            msg.setInformativeText("Connect to neaSNOM first! Click OK to connect!")
            button = msg.exec()
            if button == QMessageBox.Ok:
                self.connect_to_neasnom()
            else:
                pass
    
    def start_advanced_scan(self):
        # Check if positions were loaded an send them to worker thread
        if self.advanced_positions is not None:
            self.worker.advanced_positions = self.advanced_positions
        else:
            self.status_bar_update('No positions were loaded for advanced scan!')
            print("No positions loaded")
            msg = QMessageBox(self)
            msg.setWindowTitle("Missing point positions")
            msg.setText("Positions of measurement points are not defined")
            msg.setIcon(QMessageBox.Critical)
            msg.setStandardButtons(QMessageBox.Ok|QMessageBox.Cancel)
            buttonConnect = msg.button(QMessageBox.Ok)
            buttonConnect.setText('Choose file')
            msg.setInformativeText("Click 'Choose file' to browse for text file defining point positions")
            button = msg.exec()
            if button == QMessageBox.Ok:
                self.load_advanced_points()
                self.worker.advanced_positions = self.advanced_positions
            else:
                return None
        # Create map object and set up scan parameters
        self.advanced_map = advanced_scan()
        self.advanced_map.Npoints = len(self.advanced_positions)
        # Send the map object to worker object
        self.worker.advanced_scan_map = self.advanced_map
        # Check if connected to SNOM
        # if self.connected:
        #     # Pass SDK objects to worker thread
        #     self.worker.nea = self.nea
        #     self.worker.context = self.context
        #     self.worker.motors = self.motors
        #     self.worker.Vector3D = self.Vector3D
        # Emit Signal to start scan at worker thread Slot
        self.advanced_work_requested.emit()
            # self.connect_snom_button.setEnabled(False)
        # else:
        #     self.status_bar_update('Connect to neaSNOM before scanning!')
        #     print("No device was found")
        #     msg = QMessageBox(self)
        #     msg.setWindowTitle("No connection!")
        #     msg.setText("Not connected to SNOM!")
        #     msg.setIcon(QMessageBox.Critical)
        #     msg.setStandardButtons(QMessageBox.Ok|QMessageBox.Cancel)
        #     buttonConnect = msg.button(QMessageBox.Ok)
        #     buttonConnect.setText('Connect')
        #     msg.setInformativeText("Connect to neaSNOM first! Click OK to connect!")
        #     button = msg.exec()
        #     if button == QMessageBox.Ok:
        #         self.connect_to_neasnom()
        #     else:
        #         pass

    def abort_advanced(self):
        self.worker.advanced_aborted = True
        self.load_meas_adv_button.setEnabled(True)
        self.load_coords_button.setEnabled(True)
        self.show_coords_button.setEnabled(True)
        self.star_scan_adv_button.setEnabled(True)
        self.abort_adv_button.setEnabled(False)
            
    def update_scan_progress(self, v):
        self.currentZindex = v
        self.mirror_map = self.worker.scan_map
        self.set_display_data(self.mirror_map)
        self.update_image()
        # self.statusbar.showMessage(f'X = {self.worker.scan_map.X[-1]}, Y = {self.worker.scan_map.Y[-1]}, Z = {self.worker.scan_map.Z[-1]}')

    def update_advanced_scan_progress(self, v):
        self.advanced_map = self.worker.advanced_scan_map
        self.set_advanced_display(self.advanced_map)
        self.update_advanced_plot()

    def scan_complete(self):
        # Push measured map to display
        self.mirror_map = self.worker.scan_map
        self.center_pos_abs = self.mirror_map.center_point
        self.center_pos_rel = [0,0]
        self.center_marker = [{'pos': self.center_pos_rel, 'data': 1}]
        self.set_display_data(self.mirror_map)
        self.update_image()
        self.loaded_map = None
        self.connect_snom_button.setEnabled(True)
        if os.path.exists("temp.txt"):
            os.remove("temp.txt")
        if self.AutosaveCheckBox.isChecked():
            if self.mirror_map.Nz == 1:
                fname = f'{datetime.datetime.now().strftime("%Y.%m.%d-%H.%M")}_2D_Mirror_scan_{self.mirror_map.sizeX}x{self.mirror_map.sizeY}um.dat'
            else:
                fname = f'{datetime.datetime.now().strftime("%Y.%m.%d-%H.%M")}_3D_Mirror_scan_{self.mirror_map.sizeX}x{self.mirror_map.sizeY}x{self.mirror_map.sizeZ}um.dat'
            self.save_data(fname=fname)

    def advanced_scan_complete(self):
        # Push measured map to display
        self.advanced_map = self.worker.advanced_scan_map
        self.center_pos_abs = self.advanced_map.center_point
        self.center_pos_rel = [0,0]
        self.center_marker = [{'pos': self.center_pos_rel, 'data': 1}]
        self.set_advanced_display(self.advanced_map)
        self.update_advanced_plot()
        self.loaded_advanced_map = None

        # Enable buttons again
        # self.connect_snom_button_adv.setEnabled(True)
        self.load_meas_adv_button.setEnabled(True)
        self.load_coords_button.setEnabled(True)
        self.show_coords_button.setEnabled(True)
        self.star_scan_adv_button.setEnabled(True)
        self.abort_adv_button.setEnabled(False)

        if os.path.exists("temp.txt"):
            os.remove("temp.txt")
        if self.AutosaveCheckBox.isChecked():
            fname = f'{datetime.datetime.now().strftime("%Y.%m.%d-%H.%M")}_Nonuniform_Mirror_scan_{self.advanced_map.Npoints}point.dat'
            self.save_data(fname=fname)

    def advanced_scan_started(self):
        self.load_meas_adv_button.setEnabled(False)
        self.load_coords_button.setEnabled(False)
        self.show_coords_button.setEnabled(False)
        self.star_scan_adv_button.setEnabled(False)
        self.abort_adv_button.setEnabled(True)

    def status_bar_update(self, m):
        self.statusbar.showMessage(m)
        # print(m)

    def enable_move_to_point(self):
        if self.connected:
            if self.mirror_map is not None:
                self.click_move_enabled = True
                self.move_to_button.setEnabled(False)
            else:
                self.status_bar_update('Invalid map type for goto move!')
                print("Invalid map type for goto move!")
                msg = QMessageBox(self)
                msg.setWindowTitle("Invalid map")
                msg.setText("Invalid map!")
                msg.setIcon(QMessageBox.Critical)
                msg.setStandardButtons(QMessageBox.Close)
                msg.setInformativeText("Conduct a mirror scan first to be able to move to a specific position!")
                msg.exec()

    def save_data(self, fname):
        if self.mirror_map is not None:
            X = self.mirror_map.X.flatten()
            Y = self.mirror_map.Y.flatten()
            Z = self.mirror_map.Z.flatten()
            O1A = self.mirror_map.O1A.flatten()
            O2A = self.mirror_map.O2A.flatten()
            O3A = self.mirror_map.O3A.flatten()
            O4A = self.mirror_map.O4A.flatten()
            M = np.array([X,Y,Z,O1A,O2A,O3A,O4A])
            np.savetxt(fname, M.T,
                        header='\n'.join([f'SizeX = {self.mirror_map.sizeX}', f'SizeY = {self.mirror_map.sizeY}',f'SizeZ = {self.mirror_map.sizeZ}',
                        f'StepX = {self.mirror_map.step_sizeX}',f'StepY = {self.mirror_map.step_sizeY}',f'StepZ = {self.mirror_map.step_sizeZ}']))

    def link_scan_size(self):
        if self.sizes_linked:
            self.linkSizeButton.setStyleSheet(self.linkButton_default_style_sheet)
            self.sizeY_spinBox.setEnabled(True)
            self.sizeX_spinBox.valueChanged.disconnect()
            self.sizes_linked = False
        else:
            self.linkSizeButton.setStyleSheet("background-color: rgb(100, 100, 100); border-radius: 4px; border-color: rgb(150, 150, 150); border-width: 2px; border-style: inset; padding: 3px;")
            self.sizeY_spinBox.setEnabled(False)
            value = self.sizeX_spinBox.value()
            self.sizeY_spinBox.setValue(value)
            self.sizeX_spinBox.valueChanged.connect(self.update_linked_size)
            self.sizes_linked = True

    def link_scan_step_size(self):
        if self.step_sizes_linked:
            self.linkStepSizeButton.setStyleSheet(self.linkButton_default_style_sheet)
            self.stepY_spinBox.setEnabled(True)
            self.stepX_spinBox.valueChanged.disconnect()
            self.step_sizes_linked = False
        else:
            self.linkStepSizeButton.setStyleSheet("background-color: rgb(100, 100, 100); border-radius: 4px; border-color: rgb(150, 150, 150); border-width: 2px; border-style: inset; padding: 3px;")
            self.stepY_spinBox.setEnabled(False)
            value = self.stepX_spinBox.value()
            self.stepY_spinBox.setValue(value)
            self.stepX_spinBox.valueChanged.connect(self.update_linked_stepsize)
            self.step_sizes_linked = True

    def update_linked_stepsize(self):
            value = self.stepX_spinBox.value()
            self.stepY_spinBox.setValue(value)

    def update_linked_size(self):
            value = self.sizeX_spinBox.value()
            self.sizeY_spinBox.setValue(value)
        
class mirror_scan:
    def __init__(self):
        # Parameters
        self.center_point = None
        self.step_sizeX = None
        self.step_sizeY = None
        self.step_sizeZ = None

        self.sizeX = None
        self.sizeY = None
        self.sizeZ = None

        self.Nx = None
        self.Ny = None
        self.Nz = None

        # Data
        self.O1A = []
        self.O2A = []
        self.O3A = []
        self.O4A = []

        self.X = []
        self.Y = []
        self.Z = []

    def recalc_size(self):
        self.Nx = int(self.sizeX/self.step_sizeX)
        self.Ny = int(self.sizeY/self.step_sizeY)
        self.Nz = int(self.sizeZ/self.step_sizeZ)
        if self.Nz == 0:
            self.Nz = 1
        else:
            self.Nz = int(self.sizeZ/self.step_sizeZ)

    def create_array(self):
        self.O1A = np.zeros((self.Nz,self.Nx,self.Ny))
        self.O2A = np.zeros((self.Nz,self.Nx,self.Ny))
        self.O3A = np.zeros((self.Nz,self.Nx,self.Ny))
        self.O4A = np.zeros((self.Nz,self.Nx,self.Ny))

        self.X = np.zeros((self.Nz,self.Nx,self.Ny))
        self.Y = np.zeros((self.Nz,self.Nx,self.Ny))
        self.Z = np.zeros((self.Nz,self.Nx,self.Ny))

class advanced_scan:
    def __init__(self):
        # Parameters
        self.center_point = None

        self.Npoints = None

        # Data
        self.O1A = []
        self.O2A = []
        self.O3A = []
        self.O4A = []

        self.X = []
        self.Y = []
        self.Z = []

if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec())