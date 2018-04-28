# ------------------------------------------------------------------ #
# BG_simple.py
# ------------------------------------------------------------------ #
#
# A series of simple Bragg Gratings
#
# ------------------------------------------------------------------ #
# VERSION HISTORY
# 26 Apr 2018 - AMH - Initialization
# ------------------------------------------------------------------ #
#
# ------------------------------------------------------------------ #
#      To - DO
# ------------------------------------------------------------------ #
#
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
#      Initialize parent cell and set params
# ------------------------------------------------------------------ #

BGLattice = gdspy.Cell('BGLattice')

# GDS Parameters
layerNumber = 1        # silicon layer number

# Chip parameters
chipDim = 8.78 * 1e3   # Width / height of our chip

# Waveguide parameters
waveguideWidth      = 0.5  # width of all our waveguides
waveguideBendRadius = 10   # bend radius for waveguides

# Coupler parameters
couplePitch         = 127     # separation distance between couplers
couplerBufferMiddle = 50      # distance from drop port to taper waveguide
couplerBufferBottom = 100     # distance from trhough port to taper waveguide
middleTapperHeightBuffer = 20 # distance from middle taper bus to bottom taper bus

# Bragg parameters
gratingPeriod      = [ 0.316, 0.317, 0.318, 0.319]       # Starting ring radius
numGratPeriods     = 100       # Number of periods
gratingDuplicates  = 10         # Number of times to duplicate grating
dwidth             = [0.01, 0.02, 0.05]   # Separation between ring and busses


# Taper Parameters
taperBufferWidth = 0.180;   # width of the buffer after taper
taperBufferLength = 5;      # length of the buffer after taper

# Text parameters
textBuffer = 20;      # distance from text to edge of chip and taper
textSize   = 15;      # base size of text

# Vernier parameters
vernierPadding = 0    # Pad the corner from the floor plan

# ------------------------------------------------------------------ #
#      Add a floor plan
# ------------------------------------------------------------------ #

# Add floor plan
BGLattice.add(gdspy.Rectangle([0,0],[chipDim,-chipDim],layer=100))

# ------------------------------------------------------------------ #
#      Add Vernier Pattern
# ------------------------------------------------------------------ #

vernierCell = obLib.vernier(layerNumber)

# get dimensions of vernier pattern
vernierDims   = vernierCell.get_bounding_box()
vernierWidth  = abs(vernierDims[0,0] - vernierDims[1,0])
vernierHeight = abs(vernierDims[0,1] - vernierDims[1,1])

# Add vernier pattern
pos = (chipDim - (vernierPadding), -(vernierPadding))
BGLattice.add(gdspy.CellReference(vernierCell,pos))

# ------------------------------------------------------------------ #
#      Create several bragg gratings along chip
# ------------------------------------------------------------------ #

# Calculate the number of repeats we can fit on the chip
numPeriods = len(gratingPeriod)
numWidths = len(dwidth)
totalCombos = numPeriods * numWidths * couplePitch
print(chipDim/totalCombos)
numRepeats = np.floor((chipDim - vernierHeight)/totalCombos)    # Number of MZIs given pitch
print(numRepeats)
# Initialize coupling taper
couplingTaperCell = obLib.couplingTaper()

# Get dimensions of taper
taperDims   = couplingTaperCell.get_bounding_box()
taperWidth  = abs(taperDims[0,0] - taperDims[1,0])
taperHeight = abs(taperDims[0,1] - taperDims[1,1])

yOffset = vernierHeight;

# Iterate through
iterCoupler = 1;
for iterPeriods in range(0,int(numPeriods)):
    for iterWidths in range(0,int(numWidths)):
        # Initialize bragg grating cell
        braggCell = obLib.braggGrating(period = gratingPeriod[iterPeriods], NG = numGratPeriods,
                         waveguideWidth = waveguideWidth, dwidth = dwidth[iterWidths],layerNumber=1)

        # Initialize unit cell
        cellName = 'BraggUnit_period='+str(int(gratingPeriod[iterPeriods]*1e3))+'_width='+str(int(dwidth[iterWidths]*1e3))
        BraggUnitCell = gdspy.Cell(cellName)

        # Get dimensions of the bragg grating to align waveguides
        braggcellCDims   = braggCell.get_bounding_box()
        braggcellWidth   = abs(braggcellCDims[0,0] - braggcellCDims[1,0])
        braggcellHeight  = abs(braggcellCDims[0,1] - braggcellCDims[1,1])

        # iterate through number of repeats
        for ktaper in range(0,int(numRepeats)):

            L = gratingPeriod[iterPeriods] * numGratPeriods

            numCols = int(np.floor((chipDim - 2 * taperWidth)/L))

            # add bragg grating
            pos = (taperWidth + braggcellWidth/2,-(iterCoupler+ktaper)*couplePitch - yOffset)
            BraggUnitCell.add(gdspy.CellArray(braggCell,origin=pos,columns=numCols,rows=1, spacing=[L,0]))

            # add coupling taper
            posTaper = (taperWidth/2,-(iterCoupler+ktaper)*couplePitch - yOffset)
            BraggUnitCell.add(gdspy.CellReference(couplingTaperCell,posTaper))

            #add coupling buffer
            BraggUnitCell.add(gdspy.Rectangle(
                [-taperBufferLength,-(iterCoupler+ktaper)*couplePitch+taperBufferWidth/2 - yOffset],
                [0,-(iterCoupler+ktaper)*couplePitch-taperBufferWidth/2 - yOffset],
                layer=layerNumber))

            # add termination taper
            posTaper = (chipDim - taperWidth/2,
                        -(iterCoupler+ktaper)*couplePitch - yOffset)
            BraggUnitCell.add(gdspy.CellReference(couplingTaperCell,posTaper,rotation=180))

            #add termination buffer
            BraggUnitCell.add(gdspy.Rectangle(
                [chipDim,-(iterCoupler+ktaper)*couplePitch+taperBufferWidth/2 - yOffset],
                [chipDim + taperBufferLength,-(iterCoupler+ktaper)*couplePitch-taperBufferWidth/2 - yOffset],
                layer=layerNumber))

            # connect termination taper to bragg array
            BraggUnitCell.add(gdspy.Rectangle(
                [taperWidth + numCols * L, -(iterCoupler+ktaper)*couplePitch - yOffset + waveguideWidth/2],
                [chipDim - taperWidth, -(iterCoupler+ktaper)*couplePitch - yOffset - waveguideWidth/2],
                layer=layerNumber
            ))

            # add text below taper on left
            label = 'p='+str(int(gratingPeriod[iterPeriods]*1e3))+' w='+str(int(dwidth[iterWidths]*1e3))
            text = gdspy.Text(label, textSize,
                    (textBuffer, -(iterCoupler+ktaper)*couplePitch - textBuffer - textSize/2 - yOffset),layer=layerNumber)
            BraggUnitCell.add(text)

            # add text below taper on right
            text = gdspy.Text(label, textSize,
                    (chipDim - 8*textBuffer, -(iterCoupler+ktaper)*couplePitch - textBuffer - textSize/2 - yOffset),layer=layerNumber)
            BraggUnitCell.add(text)

        # consolidate cells to master cell
        BGLattice.add(gdspy.CellReference(BraggUnitCell))

        iterCoupler = iterCoupler + numRepeats

# ------------------------------------------------------------------ #
#      Center Parent Cell for fab
# ------------------------------------------------------------------ #
centeredCell = gdspy.Cell('centeredCell')

RRDims   = BGLattice.get_bounding_box()
RRWidth  = abs(RRDims[0,0] - RRDims[1,0])
RRHeight = abs(RRDims[0,1] - RRDims[1,1])

pos = (-RRWidth/2,RRHeight/2)
centeredCell.add(gdspy.CellReference(BGLattice,pos))
# ------------------------------------------------------------------ #
#      OUTPUT GDS FILE
# ------------------------------------------------------------------ #

# Output the layout to a GDSII file (default to all created cells).
# Set the units we used to micrometers and the precision to nanometers.

filename = 'BG_simple.gds'
outPath = os.path.abspath(os.path.join(d, 'GDS/'+filename))
gdspy.write_gds(outPath, unit=1.0e-6, precision=1.0e-9)
