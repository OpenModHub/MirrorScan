# MirrorScan [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]
Python application to realize focus spot scanning with neaSNOM microscopes.
It works only with neaSNOM devices equipped with position sensors on the mirror motors.

## Configuration

Before starting the application for the first time please enter the appropriate information in the `config.yaml` file.

Snapshot of the application window:
![app_screenshot](/images/app_screenshot.png)

## Working principle
The application records the demodulated optical signals while changing the position of the parabolic mirror across the defined area.
- the initial position of the parabolic mirror will be in the middle of the scanning area
- the scanning starts from negative coordinates for all axis
  - it measures row-by-row in the X direction 
  - when Z distance is defined it measures a 2D (X-Y) map at the defined Z positions

Schematics of the scanning directions:

![scanning](/images/scanning_schematics.png)

## Functionality
You can use the application to:
- open and display previously saved mirror scan maps
  - in the current version: :warning: you have to change Size X,Y,Z and step sizes according to the loaded measurement (will fixed soon)
- measure new mirror scans
  - after each scan, the software autosaves the resulting map in a text file containing:
    - X, Y, Z coordinates
    - O1A, O2A, O3A, O4A optical signal maps

### Before scan:
1. Make sure that the detector is cooled down and is in the right position
2. Make sure your laser (or other light source) is turned on and the focus is supposedly in nearby the tip
3. Approach to contact before starting the scan
4. You can only start scanning if you are connected to neaServer
  - use the Connect button to do so

:bulb: It is always a good practice to save the position of the mirror in neaSCAN before starting a scan

:bulb: Tipp: you can monitor the relative position of the mirror also in neaSCAN

### After scan:
While you are connected, you can move the mirror position to the desired position of the scanned area:
1. Click the `Move to` button
2. Click on a chosen position in the image
   - a small marker will move and show the new location
3. If you start a new mirror scan, this new location will be the center point of the new map
  
## Software versions
The application was tested on a device with the following software version.
- neaSCAN 2.2.10875
- neaServer 2.1.11062.0

### License

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
