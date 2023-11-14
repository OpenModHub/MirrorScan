# MirrorScan
Python application to realize focus spot scanning with neaSNOM microscopes.
It works only with neaSNOM devices equipped with position sensors on the mirror motors.

## Working principle
The application records the demodulated optical signals while changing the position of the parabolic mirror accross the define area.
- the initial position of the parabolic mirror will be the middle of the scanning area
- the scanning starts from negative coordinates for all axis
  -- it measures row-by-row in X direction
  -- when Z distance is defined it measures 2D (X-Y) map at the defined Z positions
  
## SDK version
The application was tested on a device with the following software version.
- neaSCAN 2.2.10875
- neaServer 2.1.11062.0
