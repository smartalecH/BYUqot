# ------------------------------------------------------------------ #
# doubleCouple.py
# ------------------------------------------------------------------ #
#
# A mask design used to 3D print onto a SiP chip. This design has couplers on
# both sides of the chip
#
# ------------------------------------------------------------------ #
# VERSION HISTORY
# 10 Apr 2018 - AMH - Initialization
#
# ------------------------------------------------------------------ #

# ------------------------------------------------------------------ #
#      Import libraries
# ------------------------------------------------------------------ #

# Get project library path to import library files
import sys
import os
d = os.path.dirname(os.getcwd())
libPath = os.path.abspath(os.path.join(d, 'lib'))
sys.path.insert(0, libPath)

# Import all other libraries
import gdspy
import numpy as np
import objectLibrary as obLib

# ------------------------------------------------------------------ #
#      Design Constants
# ------------------------------------------------------------------ #

layerNumber       = 1

couplerSeparation = 127      # Separation between couplers
numCouplers       = 10       # Number of couplers
taperLength       = 50       # Length of taper
taperWidth        = 0.180    # Short width of taper
waveguideWidth    = 0.5      # Width of waveguid
bufferLength      = 10       # Length of the buffer used to connect taper off chip

ringRadius        = 5        # Radius of ring resonator
ringWL            = 0.5      # Width of ring waveguides
ringBW            = 0.5      # Width of bus waveguides
racetrackLength   = 10       # Length of racetrack connecting waveguide rings
ringGap           = 0.2      # Gap between ring and bus

# ------------------------------------------------------------------ #
#      Ring resonator
# ------------------------------------------------------------------ #

def ringResonator():
    # Initialize cell
    ringCell = gdspy.Cell('Ring')

    # Calculate center ring coordinates
    ringCenterX = racetrackLength/2;
    ringCenterY = -(ringWL/2 + ringGap + waveguideWidth)/2;

    # create left ring
    ringRight = gdspy.Round((ringCenterX,
            ringCenterY), ringRadius+ringWL/2, inner_radius=ringRadius-ringWL/2,
            initial_angle=np.pi/2,final_angle=-np.pi/2,layer=layerNumber)

    # create right ring
    ringLeft = gdspy.Round((-ringCenterX,
            ringCenterY), ringRadius+ringWL/2, inner_radius=ringRadius-ringWL/2,
            initial_angle=np.pi/2,final_angle=3*np.pi/2,layer=layerNumber)

    # create top race topRaceTrack
    topRaceTrack = gdspy.Rectangle(
        [-racetrackLength/2,ringCenterY+ringRadius-ringWL/2],
        [racetrackLength/2,ringCenterY+ringRadius+ringWL/2],layer=layerNumber)

    # create bottom racetrack
    bottomRaceTrack = gdspy.Rectangle(
        [-racetrackLength/2,ringCenterY-ringRadius-ringWL/2],
        [racetrackLength/2,ringCenterY-ringRadius+ringWL/2],layer=layerNumber)

    # create top waveguide bus
    busLength = 2*ringRadius + racetrackLength
    busCenterY = ringCenterY + ringRadius + waveguideWidth/2 + ringGap + waveguideWidth
    topBus = gdspy.Rectangle([-busLength/2,busCenterY-waveguideWidth/2],
            [busLength/2,busCenterY+waveguideWidth/2],layer=layerNumber)

    # Add all geometries
    ringCell.add(topBus)
    ringCell.add(topRaceTrack)
    ringCell.add(bottomRaceTrack)
    ringCell.add(ringRight)
    ringCell.add(ringLeft)

    # return Cells
    return ringCell

# ------------------------------------------------------------------ #
#      Taper
# ------------------------------------------------------------------ #
def taper():
    # Initialize cells
    taperCell = gdspy.Cell('Taper')

    # Create geometry
    taperPts = [(-taperLength/2,-taperWidth/2),(-taperLength/2,taperWidth/2),(taperLength/2,waveguideWidth/2),(taperLength/2,-waveguideWidth/2)]
    taperPoly = gdspy.Polygon(taperPts,layer=layerNumber)

    # Add geometry
    taperCell.add(taperPoly)

    # Return cell
    return taperCell

# ------------------------------------------------------------------ #
#      Total Chip
# ------------------------------------------------------------------ #
def chip():
    # Initialize cell
    chipCell = gdspy.cell('doubleCouple')

    # Iterate through all couplers

# ------------------------------------------------------------------ #
#      OUTPUT
# ------------------------------------------------------------------ #


# Output the layout to a GDSII file (default to all created cells).
# Set the units we used to micrometers and the precision to nanometers.
filename = 'doubleCouple.gds'
outPath = os.path.abspath(os.path.join(d, 'GDS/'+filename))
gdspy.write_gds(outPath, unit=1.0e-6, precision=1.0e-9)
