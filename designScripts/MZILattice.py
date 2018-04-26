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
#      Initialize parent cell
# ------------------------------------------------------------------ #

MZILatticeCell = gdspy.Cell('MZILattice')

# ------------------------------------------------------------------ #
#      Add a floor plan
# ------------------------------------------------------------------ #

# Width / height of our chip
chipDim = 8.78 * 1e3

# width of all our waveguides
waveguideWidth = 0.5;

deltaL = 0
layerNumber = 1
couplerBufferMiddle = 50
couplerBufferBottom = 100
TMradius = 10

repeatDelta = 3;
Lincrement  = 100;

# Text parameters
textBuffer = 20;      # distance from text to edge of chip and taper
textSize   = 15;      # base size of text

# Add floor plan
MZILatticeCell.add(gdspy.Rectangle([0,0],[chipDim,-chipDim],layer=100))

# ------------------------------------------------------------------ #
#      Create several MZI's along chip
# ------------------------------------------------------------------ #

couplePitch = 127     # distance between couplers
numMZI = np.floor(chipDim/couplePitch)    # Number of MZIs given pitch

# Initialize two kinds of MZIs and taper
#MZIcellY = obLib.MZI(deltaL = deltaL, bendRadius = TMradius, Lref = 100, gapLength = 50, waveguideWidth = 0.5,coupleType="Y");

couplingTaperCell = obLib.couplingTaper()

# Get dimensions of taper
taperDims   = couplingTaperCell.get_bounding_box()
taperWidth  = abs(taperDims[0,0] - taperDims[1,0])
taperHeight = abs(taperDims[0,1] - taperDims[1,1])

# Iterate through MZIs and draw system
incrementOffset = 0;
for k in range(1,int(numMZI)):
    if k % 3 == 0:
        if incrementOffset % repeatDelta == 0:
            deltaL = deltaL + Lincrement;
            MZIcellC = obLib.MZI(deltaL = deltaL, bendRadius = TMradius,
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
        l1path.fillet(radius=TMradius)
        MZIUnitCell.add(l1path)

        # connect bottom taper to coupler
        length = [k*couplePitch - taperWidth + MZIcellWidth/2 + couplerBufferBottom,
                  2*couplePitch + 2*yOffset - waveguideWidth, couplerBufferBottom + 1]
        turn = [1,1]
        l2path = gdspy.L1Path(initial_point=(taperWidth, -(k+2)*couplePitch),
                              direction='+x', width=waveguideWidth, length=length,
                              turn=turn, layer=layerNumber)
        l2path.fillet(radius=TMradius)
        MZIUnitCell.add(l2path)

        # consolidate cells to master cell
        MZILatticeCell.add(gdspy.CellReference(MZIUnitCell))

        incrementOffset = incrementOffset + 1

# ------------------------------------------------------------------ #
#      Vernier Pattern
# ------------------------------------------------------------------ #

vernierCell = obLib.vernier(1)

# get dimensions of vernier pattern
vernierDims   = vernierCell.get_bounding_box()
vernierWidth  = abs(vernierDims[0,0] - vernierDims[1,0])
vernierHeight = abs(vernierDims[0,1] - vernierDims[1,1])

# Padd the corner from the floor plan
padding = 30
# Add vernier pattern
pos = (chipDim - (padding), -(padding))
MZILatticeCell.add(gdspy.CellReference(vernierCell,pos))

# ------------------------------------------------------------------ #
#      OUTPUT
# ------------------------------------------------------------------ #

# Output the layout to a GDSII file (default to all created cells).
# Set the units we used to micrometers and the precision to nanometers.

filename = 'MZI.gds'
outPath = os.path.abspath(os.path.join(d, 'GDS/'+filename))
gdspy.write_gds(outPath, unit=1.0e-6, precision=1.0e-9)
