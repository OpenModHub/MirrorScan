# MirrorScan
Python application to realize focus spot scanning with neaSNOM microscopes.
It works only with neaSNOM devices equipped with position sensors on the mirror motors.

## Working principle
The application records the demodulated optical signals while changing the position of the parabolic mirror accross the define area.
- the initial position of the parabolic mirror will be the middle of the scanning area
- the scanning starts from negative coordinates for all axis
  - it measures row-by-row in X direction
   when Z distance is defined it measures 2D (X-Y) map at the defined Z positions

## Functionality
Use can use the application to:
- open and display previously saved mirror scan maps
- measure new mirror scans
  - after each scan the software autosaves the resulting map in a text files containing:
    - X, Y, Z coordinates
    - O1A, O2A, O3A, O4A optical signal maps

### Before mirror scan:
1. Make sure that the detector is cooled down and is in the right position
2. Make sure your laser (or other light source) is turned on and the focus is supposedly in nearby the tip
3. Approach to contact before starting the scan
4. You can only start scanning if you are connected to neaServer
  - use the Connect button to do so

:bulb: It is always a good practice to save the position of the mirror in neaSCAN before starting a scan

:bulb: Tipp: you can monitor the relative position of the mirror also in neaSCAN

### After scan:
While your are connected, you can move the mirror position to a desired position of the scanned area:
1. Click "move to" button
2. Click in the image
   - a small marker will move and show the new location
3. If you start a new mirror scan, this new location will be the center point of the new map
  
## SDK version
The application was tested on a device with the following software version.
- neaSCAN 2.2.10875
- neaServer 2.1.11062.0
