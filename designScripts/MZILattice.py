# ------------------------------------------------------------------ #
# MZILattice.py
# ------------------------------------------------------------------ #
#
# A series of MZI's for 3D printing testing purposes
#
# ------------------------------------------------------------------ #
# VERSION HISTORY
# 10 Apr 2018 - AMH - Initialization
# 24 Apr 2018 - AMH - Added vernier pattern, directional couplers
# 25 Apr 2018 - AMH - Varied the deltaL for each MZI
# ------------------------------------------------------------------ #
#
# ------------------------------------------------------------------ #
#      To - DO
# ------------------------------------------------------------------ #
# Code cleanup
# Add bosch etch buffer off of floor plan
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
#      Set params
# ------------------------------------------------------------------ #

# Initalize parent cell
MZILatticeCell = gdspy.Cell('MZILattice')

# Chip params
layerNumber = 1              # silicon layer
chipDim     = 8.78 * 1e3     # Width / height of our chip
couplePitch = 127            # distance between couplers

# Waveguide parameters
waveguideWidth      = 0.5    # width of all our waveguides
waveguideBendRadius = 10     # bend radius of all waveguides in microns

# MZI parameters
deltaL      = 0      # Initial change in path length
repeatDelta = 3      # Number of times to repeat the same MZI
Lincrement  = 100    # Step wise increase in MZI

# Coupler parameters
couplerBufferMiddle = 50   # Separation between bottom coupler and waveguide bend
couplerBufferBottom = 100  # Separation between top coupler and waveguide bend

# Taper Parameters
taperBufferWidth = 0.180;   # width of the buffer after taper
taperBufferLength = 5;      # length of the buffer after taper

# Text parameters
textBuffer = 20;      # distance from text to edge of chip and taper
textSize   = 15;      # base size of text

# Vernier parameters
vernierPadding = 0   # Padd the corner from the floor plan

# ------------------------------------------------------------------ #
#      Add a floor plan
# ------------------------------------------------------------------ #

MZILatticeCell.add(gdspy.Rectangle([0,0],[chipDim,-chipDim],layer=100))

# ------------------------------------------------------------------ #
#      Create several MZI's along chip
# ------------------------------------------------------------------ #

# Number of MZIs given pitch
numMZI = np.floor(chipDim/couplePitch)

# Initialize taper cell
couplingTaperCell = obLib.couplingTaper()

# Get dimensions of taper
taperDims   = couplingTaperCell.get_bounding_box()
taperWidth  = abs(taperDims[0,0] - taperDims[1,0])
taperHeight = abs(taperDims[0,1] - taperDims[1,1])

# Iterate through MZIs and draw system
incrementOffset = 0;
for k in range(1,int(numMZI)):

    # Iterate through all 3 tapers for each MZI
    if k % 3 == 0:

        # Only create new MZI cell if it's a new path length difference
        if incrementOffset % repeatDelta == 0:
            deltaL = deltaL + Lincrement;
            MZIcellC = obLib.MZI(deltaL = deltaL, bendRadius = waveguideBendRadius,
            Lref = 250, gapLength = 50, waveguideWidth = 0.5,coupleType="C");

        # Initialize unit cell
        MZIUnitCell = gdspy.Cell('MZIUnit_'+str(k))

        # Get dimensions of the directional coupler to compensate for height offset
        MZIcellCDims   = MZIcellC.get_bounding_box()
        MZIcellWidth  = abs(MZIcellCDims[0,0] - MZIcellCDims[1,0])
        yOffset = abs(MZIcellCDims[0,1])

        # add coupler
        pos = (k*couplePitch,-k*couplePitch + yOffset - waveguideWidth/2)
        MZIUnitCell.add(gdspy.CellReference(MZIcellC,pos))

        # add termination taper to prevent reflections
        posTaper = (k*couplePitch - MZIcellWidth/2 - taperWidth/2,-k*couplePitch + 2*yOffset - waveguideWidth)
        MZIUnitCell.add(gdspy.CellReference(couplingTaperCell,posTaper))

        # add coupling tapers
        for ktaper in range(0,3):
            posTaper = (taperWidth/2,-(k+ktaper)*couplePitch)
            MZIUnitCell.add(gdspy.CellReference(couplingTaperCell,posTaper))

            MZIUnitCell.add(gdspy.Rectangle(
                [-taperBufferLength,-(k+ktaper)*couplePitch+taperBufferWidth/2],
                [0,-(k+ktaper)*couplePitch-taperBufferWidth/2],
                layer=layerNumber))

        # add text below coupling first taper
        text = gdspy.Text(str(int(deltaL)), textSize,
                (textBuffer, -k*couplePitch - textBuffer - textSize/2),layer=layerNumber)
        MZIUnitCell.add(text)

        # connect top taper to coupler
        MZIUnitCell.add(gdspy.Rectangle(
            [taperWidth,-k*couplePitch + waveguideWidth/2],
            [k*couplePitch - MZIcellWidth/2,-k*couplePitch - waveguideWidth/2],layer=layerNumber))

        # connect middle taper to coupler
        length = [k*couplePitch - taperWidth + MZIcellWidth/2 + couplerBufferMiddle,
                  couplePitch, couplerBufferMiddle + 1]
        turn = [1,1]
        l1path = gdspy.L1Path(initial_point=(taperWidth, -(k+1)*couplePitch),
                              direction='+x', width=waveguideWidth, length=length,
                              turn=turn, layer=layerNumber)
        l1path.fillet(radius=waveguideBendRadius)
        MZIUnitCell.add(l1path)

        # connect bottom taper to coupler
        length = [k*couplePitch - taperWidth + MZIcellWidth/2 + couplerBufferBottom,
                  2*couplePitch + 2*yOffset - waveguideWidth, couplerBufferBottom + 1]
        turn = [1,1]
        l2path = gdspy.L1Path(initial_point=(taperWidth, -(k+2)*couplePitch),
                              direction='+x', width=waveguideWidth, length=length,
                              turn=turn, layer=layerNumber)
        l2path.fillet(radius=waveguideBendRadius)
        MZIUnitCell.add(l2path)

        # consolidate cells to master cell
        MZILatticeCell.add(gdspy.CellReference(MZIUnitCell))

        incrementOffset = incrementOffset + 1

# ------------------------------------------------------------------ #
#      Vernier Pattern
# ------------------------------------------------------------------ #

vernierCell = obLib.vernier(layerNumber)

# get dimensions of vernier pattern
vernierDims   = vernierCell.get_bounding_box()
vernierWidth  = abs(vernierDims[0,0] - vernierDims[1,0])
vernierHeight = abs(vernierDims[0,1] - vernierDims[1,1])

# Add vernier pattern
pos = (chipDim - (vernierPadding), -(vernierPadding))
MZILatticeCell.add(gdspy.CellReference(vernierCell,pos))

# ------------------------------------------------------------------ #
#      Center Parent Cell for fab
# ------------------------------------------------------------------ #
centeredCell = gdspy.Cell('centeredCell')

MZIDims   = MZILatticeCell.get_bounding_box()
MZIWidth  = abs(MZIDims[0,0] - MZIDims[1,0])
MZIHeight = abs(MZIDims[0,1] - MZIDims[1,1])

pos = (-MZIWidth/2,MZIHeight/2)
centeredCell.add(gdspy.CellReference(MZILatticeCell,pos))


# ------------------------------------------------------------------ #
#      OUTPUT
# ------------------------------------------------------------------ #

# Output the layout to a GDSII file (default to all created cells).
# Set the units we used to micrometers and the precision to nanometers.

filename = 'MZI.gds'
outPath = os.path.abspath(os.path.join(d, 'GDS/'+filename))
gdspy.write_gds(outPath, unit=1.0e-6, precision=1.0e-9)
