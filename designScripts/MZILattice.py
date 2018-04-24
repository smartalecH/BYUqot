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

layerNumber = 1

# Add floor plan
MZILatticeCell.add(gdspy.Rectangle([0,0],[chipDim,-chipDim],layer=100))

# ------------------------------------------------------------------ #
#      Create several MZI's along chip
# ------------------------------------------------------------------ #

couplePitch = 250     # distance between couplers
numMZI = np.floor(chipDim/couplePitch)    # Number of MZIs given pitch

# Initialize two kinds of MZIs and taper
MZIcellY = obLib.MZI(deltaL = 200, bendRadius = 10, Lref = 100, gapLength = 50, waveguideWidth = 0.5,coupleType="Y");
MZIcellC = obLib.MZI(deltaL = 200, bendRadius = 10,  Lref = 100, gapLength = 50, waveguideWidth = 0.5,coupleType="C");
couplingTaperCell = obLib.couplingTaper()

# Get dimensions of the directional coupler to compensate for height offset
MZIcellCDims   = MZIcellC.get_bounding_box()
MZIcellWidth  = abs(MZIcellCDims[0,0] - MZIcellCDims[1,0])
yOffset = abs(MZIcellCDims[0,1])


# Get dimensions of taper
taperDims   = couplingTaperCell.get_bounding_box()
taperWidth  = abs(taperDims[0,0] - taperDims[1,0])
taperHeight = abs(taperDims[0,1] - taperDims[1,1])

# Iterate through MZIs and draw system
for k in range(1,int(numMZI)):
    if k % 3 == 0:
        # Initialize unit cell
        MZIUnitCell = gdspy.Cell('MZIUnit_'+str(k))

        # add coupler
        pos = (k*couplePitch,-k*couplePitch + yOffset - waveguideWidth/2)
        MZIUnitCell.add(gdspy.CellReference(MZIcellC,pos))

        # add tapers
        for ktaper in range(0,3):
            posTaper = (taperWidth/2,-(k+ktaper)*couplePitch)
            MZIUnitCell.add(gdspy.CellReference(couplingTaperCell,posTaper))

        # connect top taper to coupler
        MZIUnitCell.add(gdspy.Rectangle(
            [taperWidth,-k*couplePitch + waveguideWidth/2],
            [k*couplePitch - MZIcellWidth/2,-k*couplePitch - waveguideWidth/2],layer=layerNumber))

        # connect middle taper to coupler

        # connect bottom taper to coupler

        # consolidate cells to master cell
        MZILatticeCell.add(gdspy.CellReference(MZIUnitCell))

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
