# ------------------------------------------------------------------ #
# vernierMask.py
# ------------------------------------------------------------------ #
#
# A mask design used to align the 3D printer to a silicon photonic chip
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

# Cell parameters
layerNumber = 1

# Vernier mask design parameters (all values in microns)
numFingers        = 10       # Number of fingers to have on top and bottom
fingerWidth       = 30       # Width of each finger
fingerSpacing     = 40       # Spacing between fingers
longFingerLength  = 200;     # Length of the long, middle finger
shortFingerLength = 150;     # Length of the short, outer fingers
baseThickness     = 76;      # Thickness of edge border of design

separationDistance = 380     # distance from edge of pattern to origin

buffer        = 50                         # Kerf width of blade
innerBoxWidth = 8.78e3                     # Actual dimensions of chip
outerBoxWidth = innerBoxWidth + buffer     # Buffered chip size

numCells = 12             # number of repeated cells in each dimension

# Now create a series of functions that return a cell. We'll leverage the recursive
# nature of GDS files to keep things simple.

# ------------------------------------------------------------------ #
#      Create single Vernier pattern
# ------------------------------------------------------------------ #

def vernier():
    # Intialize cell
    vernierCell = gdspy.Cell('vernier')

    # Calculate properties
    vernierWidth  = (longFingerLength + baseThickness)
    vernierHeight = (2*numFingers + 1) * fingerWidth + fingerSpacing * 2 * numFingers
    xCenter       = -(vernierWidth/2 - baseThickness)

    # First, place the long finger in the middle
    vernierCell.add(gdspy.Rectangle([xCenter, -fingerWidth/2],
        [xCenter+longFingerLength, fingerWidth/2],layer=layerNumber))

    # Next, iterate through and place the other fingers
    for k in range(1,numFingers+1):
        # Add top fingers
        vernierCell.add(gdspy.Rectangle(
            [xCenter, fingerWidth/2 + fingerSpacing*k + fingerWidth*(k-1)],
            [xCenter+shortFingerLength,fingerWidth/2 + fingerSpacing*k +fingerWidth*k],
            layer=layerNumber))

        # Add bottom fingers
        vernierCell.add(gdspy.Rectangle(
            [xCenter, -(fingerWidth/2 + fingerSpacing*k + fingerWidth*(k-1))],
            [xCenter+shortFingerLength,-(fingerWidth/2 + fingerSpacing*k +fingerWidth*k)],
            layer=layerNumber))

    # Finally, add the edge
    baseHeight = (2*numFingers + 1) * fingerWidth + fingerSpacing * 2 * numFingers
    vernierCell.add(gdspy.Rectangle([-vernierWidth/2, baseHeight/2],
        [xCenter, -baseHeight/2],layer=layerNumber))

    # Now let's flatten
    vernierCell.flatten()

    # Return the cell
    return vernierCell

# ------------------------------------------------------------------ #
#      Create 2D Vernier pattern from single pattern
# ------------------------------------------------------------------ #

def vernier2D():
    # Intialize 2D cell
    vernier2DCell = gdspy.Cell('vernier2D')

    # Initialize 1D cell
    vernierCell = vernier()

    # Get vernier dimensions
    vernierDims   = vernierCell.get_bounding_box()
    vernierWidth  = abs(vernierDims[0,0] - vernierDims[1,0])
    vernierHeight = abs(vernierDims[0,1] - vernierDims[1,1])

    # Place one Vernier pattern in the x direction
    xCell = gdspy.CellReference(vernierCell,rotation=-90)
    xCell.translate(-(vernierHeight/2 + separationDistance),-vernierWidth/2)
    vernier2DCell.add(xCell)

    # Place another Vernier pattern in the y direction
    yCell = gdspy.CellReference(vernierCell,rotation=180)
    yCell.translate(-vernierWidth/2,-(vernierHeight/2 + separationDistance))
    vernier2DCell.add(yCell)

    # Return final cell
    return vernier2DCell

# ------------------------------------------------------------------ #
#      Create Box outline
# ------------------------------------------------------------------ #

def boxOutline():
    # initialize cell
    outlineCell = gdspy.Cell('outline')

    # define an outer box
    outerBox = gdspy.Rectangle([-outerBoxWidth/2,-outerBoxWidth/2],
        [outerBoxWidth/2,outerBoxWidth/2],layer=layerNumber)

    # define an inner box
    innerBox = gdspy.Rectangle([-innerBoxWidth/2,-innerBoxWidth/2],
        [innerBoxWidth/2,innerBoxWidth/2],layer=layerNumber)

    # now subtract the two
    outline = gdspy.fast_boolean(outerBox,innerBox,'xor',layer=layerNumber)

    # update the cell
    outlineCell.add(outline)

    # return the cell
    return outlineCell

# ------------------------------------------------------------------ #
#      Create Single Chip
# ------------------------------------------------------------------ #

def vernierChip():
    # Initialize cells
    vernierChipCell = gdspy.Cell('vernierChip')
    vernier2DCell   = vernier2D()
    boxOutlineCell  = boxOutline()

    # Add border first
    vernierChipCell.add(gdspy.CellReference(boxOutlineCell,(0,0)))
    chipDims   = vernierChipCell.get_bounding_box()
    chipWidth  = abs(chipDims[0,0] - chipDims[1,0])

    # Now iterate through placing corners
    thetaPos = [45, 135, -135, -45]
    thetaRot = [0, 90, 180, -90]
    for k in range(0,4):
        xPos = np.sign(np.cos(np.deg2rad(thetaPos[k]))) * (chipWidth/2 - buffer/2)
        yPos = np.sign(np.sin(np.deg2rad(thetaPos[k]))) * (chipWidth/2 - buffer/2)
        vernierChipCell.add(gdspy.CellReference(vernier2DCell,(xPos,yPos),rotation=thetaRot[k]))

    # return cell
    return vernierChipCell

# ------------------------------------------------------------------ #
#      Tapeout entire wafer
# ------------------------------------------------------------------ #

def vernierMask():

    # Initialize cells
    vernierMaskCell = gdspy.Cell('vernierMask')
    vernierChipCell = vernierChip()

    # Get chip dimensions
    chipDims   = vernierChipCell.get_bounding_box()
    chipWidth  = abs(chipDims[0,0] - chipDims[1,0])

    # Get mask center
    center = (numCells * chipWidth) / 2

    # Let's make an array
    vernierMaskCell.add(gdspy.CellArray(
        vernierChipCell, numCells, numCells, (chipWidth, chipWidth), (-center, -center)
    ))

    # return final cell
    return vernierMaskCell

# ------------------------------------------------------------------ #
#      OUTPUT
# ------------------------------------------------------------------ #

vernierMask()

# Output the layout to a GDSII file (default to all created cells).
# Set the units we used to micrometers and the precision to nanometers.
filename = 'vernierMask.gds'
outPath = os.path.abspath(os.path.join(d, 'GDS/'+filename))
gdspy.write_gds(outPath, unit=1.0e-6, precision=1.0e-9)
