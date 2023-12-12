import sys
from PySide6.QtWidgets import QApplication, QFileDialog, QLabel, QVBoxLayout, QWidget, QProgressBar, QMessageBox
from PySide6.QtCore import QObject, QThread, Signal, Slot
from PySide6.QtGui import QTransform
import pyqtgraph as pg
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
    print("nea_tools module not found, working in reader mode")
    offline_mode = True

current_folder = os.getcwd()
ui_file = os.path.join(current_folder,'mirrorApp.ui')

uiclass, baseclass = pg.Qt.loadUiType(ui_file)

######## QT WORKING THREAD CLASS ############
class Worker(QObject):
    progress = Signal(int)
    completed = Signal()
    started = Signal()
    status_update = Signal(str)
    def __init__(self):
        super().__init__()
        self.scan_map = []
        self.nea = None
        self.context = None
        self.motors = None
        self.Vector3D = None

    @Slot()
    def do_scan_test(self):
        self.started.emit()
        # Calculate mirror coordinates for movement
        xs = np.linspace(-self.scan_map.sizeX/2,self.scan_map.sizeX/2,self.scan_map.Nx)
        ys = np.linspace(-self.scan_map.sizeY/2,self.scan_map.sizeY/2,self.scan_map.Ny)
        if self.scan_map.Nz == 1:
            zs = np.array([0])
        else:
            zs = np.linspace(-self.scan_map.sizeZ/2,self.scan_map.sizeZ/2,self.scan_map.Nz)

        print(f"Nz in the working thread scan: {self.scan_map.Nz}")

        # SCANNING LOOP
        counter = 0
        for idz, z in enumerate(zs):
            for idy, y in enumerate(ys):
                for idx, x in enumerate(xs):
                    counter += 1
                    self.scan_map.O1A[idz,idx,idy] = np.random.rand()
                    self.scan_map.O2A[idz,idx,idy] = np.random.rand()
                    self.scan_map.O3A[idz,idx,idy] = np.random.rand()
                    self.scan_map.O4A[idz,idx,idy] = np.random.rand()
                    self.scan_map.X[idz,idx,idy] = x
                    self.scan_map.Y[idz,idx,idy] = y
                    self.scan_map.Z[idz,idx,idy] = z
                    sleep(0.1)
                    self.progress.emit(counter)
        self.completed.emit()

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
                    self.progress.emit(counter)
                    self.status_update.emit(f'X: {newx}, Y: {newy}, Z: {newz} Remaining time: {datetime.timedelta(seconds=remtime)}')
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

######## MAIN APPLICATION WINDOW CLASS ############
class MainWindow(uiclass, baseclass):
    work_requested = Signal()
    pg.setConfigOptions(imageAxisOrder='row-major')

    def __init__(self):
        super().__init__()

        # Load UI
        self.setupUi(self)

        # Stylize
        self.setWindowTitle('Focus scanner application')
        # self.setStyleSheet("background-color: white;")
        self.channel_comboBox.addItems(['O1A', 'O2A', 'O3A','O4A'])
        self.channel_comboBox.setCurrentText('O3A')

        # Linking button label correction
        txt = "\U0001F517"
        self.linkSizeButton.setText(txt)
        self.linkButton_default_style_sheet = self.linkSizeButton.styleSheet()

        # Create the worker thread
        self.worker = Worker()
        self.worker_thread = QThread()
        self.worker.progress.connect(self.update_scan_progress)
        self.worker.completed.connect(self.scan_complete)
        self.worker.status_update.connect(self.status_bar_update)
        self.work_requested.connect(self.worker.do_scan)
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
        self.plot_area.setBackground('w')
        self.plot_area.setMouseTracking(True)                                                       # For cursor tracking
        self.imItem.hoverEvent = self.imageHoverEvent                                               # Attach event
        self.imItem.mouseClickEvent = self.imageClickEvent                                          # Attach event

        # Connect button signals
        self.choose_file_button.clicked.connect(self.choose_file)
        self.datascroll_spinBox.valueChanged.connect(lambda: self.data_scroll())
        self.channel_comboBox.currentIndexChanged.connect(self.channel_change)
        self.linkSizeButton.clicked.connect(self.link_scan_size)

        if not offline_mode:
            self.scan_button.clicked.connect(self.start_scan)
            self.connect_snom_button.clicked.connect(self.connect_to_neasnom)
            self.move_to_button.clicked.connect(self.enable_move_to_point)
            self.save_button.clicked.connect(self.save_data)

        # Other attributes and flags
        self.connected = False
        self.click_move_enabled = False
        self.mirror_map = None
        self.loaded_map = None
        self.sizes_linked = False

        if offline_mode:
            self.statusbar.showMessage(u"\u26A0 nea_tools module not found, running in display-only mode.")
            self.connect_snom_button.setEnabled(False)
            self.move_to_button.setEnabled(False)
            self.scan_button.setEnabled(False)
            self.save_button.setEnabled(False)

    def connect_to_neasnom(self):
        if "nea_tools" not in sys.modules:
            return

        path_to_dll = r"\\nea-server\updates\Application Files\neaSCAN_2_1_10694_0"
        fingerprint = 'af3b0d0f-cdbb-4555-9bdb-6fe200b64b51'
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
                print("Disable the WiFi connection!")
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
                self.statusbar.showMessage("Connected to SNOM")
                self.connect_snom_button.setText("Disconnect from neaSNOM")

            self.context = context
            self.nea = nea
            self.motors = motors
            return context, nea

    def choose_file(self):
        fname = QFileDialog.getOpenFileName(self, "Choose file","","Datatext files (*.txt *.dat)")
        self.file_name = fname[0]
        try:
            self.load_data()
            self.statusbar.showMessage(f'self.channel_comboBox.currentText() is loaded from {fname}')
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
            number = int(header_line[idx+2:])
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

    def set_display_data(self,map):
        self.meas_data = np.array(getattr(map,self.channel_comboBox.currentText()))
        # self.meas_data = self.meas_data.reshape((map.Nz,map.Nx,map.Ny))
        if map.Nz == 1:
            self.data_to_plot = self.meas_data[0,:,:]
            self.datascroll_spinBox.setValue(0)
            self.datascroll_spinBox.setRange(0,0)
            self.Zplane_label.setText(f"This is only a 2D scan")
            self.cbar.setLevels(values = (np.min(self.data_to_plot),np.max(self.data_to_plot)))
        else:
            self.zeroZ_data = self.meas_data[int(map.Nz/2),:,:]
            self.data_to_plot = self.zeroZ_data
            self.datascroll_spinBox.setValue(int(map.Nz/2))
            self.datascroll_spinBox.setRange(0, map.Nz-1)
            self.Zplane_label.setText(f"Displayed Z plane: {self.Zaxis[int(map.Nz/2)]} nm")
            self.cbar.setLevels(values = (np.min(self.zeroZ_data),np.max(self.zeroZ_data)))

        #Set up the axis values by transforming the image
        tr = QTransform()                                                               # prepare ImageItem transformation:
        tr.translate(-map.sizeX/2,-map.sizeY/2)                                         # move 3x3 image to locate center at axis origin
        tr.scale(map.sizeX/map.Nx, map.sizeY/map.Ny)                                    # scale horizontal and vertical axes
        self.imItem.setTransform(tr)
        self.plot_area.getAxis('bottom').setLabel('X position (nm)')
        self.plot_area.getAxis('left').setLabel('Y position (nm)')
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
        if (self.loaded_map == None) | (self.mirror_map == None):
            pass
        else:
            index = self.datascroll_spinBox.value()
            self.data_to_plot = self.meas_data[index,:,:]
            self.update_image()
            self.Zplane_label.setText(f"Displayed Z plane: {self.Zaxis[index]} nm")
    
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
            p.go_relative(dx,dy,0)
            p.await_movement()
            self.status_bar_update(f"Relative move to {[x,y]}")
            # Check new position
            self.context.Microscope.RefreshActiveMotorPositionXyzAsync().Wait()
            self.center_pos_abs = p.absolute_position
            self.status_bar_update(f"Absolute AFTER center move: {self.center_pos_abs}")
            # Replace center marker
            realdx = self.center_pos_abs[0] - pos_before[0]
            realdy = self.center_pos_abs[1] - pos_before[1]
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
            button = msg.exec_()
            if button == QMessageBox.Ok:
                self.connect_to_neasnom()
            else:
                pass
            
    def update_scan_progress(self, v):
        self.mirror_map = self.worker.scan_map
        self.set_display_data(self.mirror_map)
        self.update_image()
        # self.statusbar.showMessage(f'X = {self.worker.scan_map.X[-1]}, Y = {self.worker.scan_map.Y[-1]}, Z = {self.worker.scan_map.Z[-1]}')
        
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
        if self.AutosaveCheckbox.isChecked():
            self.save_data()

    def status_bar_update(self, m):
        self.statusbar.showMessage(m)
        print(m)

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
                msg.exec_()

    def save_data(self):
        if self.mirror_map is not None:
            fname = f'{datetime.datetime.now().strftime("%Y.%m.%d-%H.%M")}_2D_Mirror_scan_{self.mirror_map.sizeX}x{self.mirror_map.sizeY}_{self.mirror_map.step_sizeX}um.dat'
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
            self.stepY_spinBox.setEnabled(True)
            self.stepX_spinBox.valueChanged.disconnect(self.update_linked_stepsize)
            self.sizes_linked = False
        else:
            self.linkSizeButton.setStyleSheet("background-color: rgb(100, 100, 100); border-radius: 4px; border-color: rgb(150, 150, 150); border-width: 2px; border-style: inset; padding: 3px;")
            self.stepY_spinBox.setEnabled(False)
            value = self.stepX_spinBox.value()
            self.stepY_spinBox.setValue(value)
            self.stepX_spinBox.valueChanged.connect(self.update_linked_stepsize)
            self.sizes_linked = True

    def update_linked_stepsize(self):
            value = self.stepX_spinBox.value()
            self.stepY_spinBox.setValue(value)
        
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

if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec())