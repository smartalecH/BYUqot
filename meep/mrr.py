# ------------------------------------------------------------------ #
#       mrr.py
# ------------------------------------------------------------------ #
#
# Extract all the relevant parameters of a microring resonator:
#   Bus effective index
#   Ring effective index
#   Self-coupling coefficient
#   Cross-coupling coefficient
#   Group index
# ------------------------------------------------------------------ #
#       VERSION HISTORY
# ------------------------------------------------------------------ #
# 16 May 2018 - AMH - Initialization
#
# ------------------------------------------------------------------ #
#       To - DO
# ------------------------------------------------------------------ #
#
# ------------------------------------------------------------------ #
# Import packages
# ------------------------------------------------------------------ #

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------ #
# Define variables
# ------------------------------------------------------------------ #

Si   = mp.Medium(epsilon=12)
SiO2 = mp.Medium(epsilon=3)

# All units in microns
ringRadius        = 5      # radius of the center of the ring waveguide
ringWaveguideWith = 0.5     # width of the ring waveguide
gapLength         = 0.2     # distance between ring edge and bus
topBusWidth       = 0.5     # width of top bus
bottomBusWidth    = 0.5     # width of bottom bus
topBus            = True    # Boolean - true if including top bus
bottomBus         = True    # Boolean - true if including bottom bus
racetrackLength   = 0       # length of the racetrack

waveguideThickness = 0   # thickness of all waveguides
sim2D = mp.inf
# ------------------------------------------------------------------ #
# Draw geometry
# ------------------------------------------------------------------ #

# preallocate geometry
geometry = [];

# draw the ring
outerRingGeom = mp.Cylinder(
                    radius = ringRadius + ringWaveguideWith / 2,
                    height = sim2D,
                    center = mp.Vector3(0,0,0),
                    material = Si
                    )
innerRingGeom = mp.Cylinder(
                    radius = ringRadius - ringWaveguideWith / 2,
                    height = sim2D,
                    center = mp.Vector3(0,0,0),
                    material = SiO2
                    )
geometry.append(outerRingGeom)
geometry.append(innerRingGeom)

# draw the top bus
if topBus:
    yPos = ringRadius + ringWaveguideWith/2 + gapLength + topBusWidth/2;
    topBusGeom = mp.Block(
       mp.Vector3(3*ringRadius,bottomBusWidth,sim2D),
       center = mp.Vector3(0,yPos,0),
       material = Si
       )
    geometry.append(topBusGeom)


# draw the bottom bus
if bottomBus:
    yPos = ringRadius + ringWaveguideWith/2 + gapLength + bottomBusWidth/2;
    bottomBusGeom = mp.Block(
       mp.Vector3(3*ringRadius,bottomBusWidth,sim2D),
       center = mp.Vector3(0,-yPos,0),
       material = Si
       )
    geometry.append(bottomBusGeom)

# ------------------------------------------------------------------ #
# Setup simulation domain
# ------------------------------------------------------------------ #

# Size of domain
bufferY = 2
cellX = 3*ringRadius
cellY = 2*(ringRadius + ringWaveguideWith/2 + gapLength + ringWaveguideWith + bufferY)
cellZ = 0
cell = mp.Vector3(cellX,cellY,cellZ)

# Sources

cY = -(ringRadius + ringWaveguideWith + gapLength)
sources = [ mp.EigenModeSource(src=mp.ContinuousSource(wavelength=1.55),
                               eig_band=1,
                               size=mp.Vector3(0,3,0),
                               center=mp.Vector3(-ringRadius,cY,0),
                               eig_parity=mp.ODD_Z
                               ) ]

# Boundary conditions
pml_layers = [mp.PML(2)]

# Resolution
resolution = 10

sim = mp.Simulation(cell_size = cell,
                    boundary_layers = pml_layers,
                    geometry = geometry,
                    sources = sources,
                    resolution = resolution,
                    default_material=SiO2,
                    force_complex_fields=True)
sim.run(until=200)
eps_data = sim.get_array(center=mp.Vector3(),size=cell,component=mp.Dielectric)
cellSize = eps_data.shape
cellSizeX = cellSize[0]
cellSizeY = cellSize[1]
#cellSizeZ = cellSize[2]
#eps_data_mod = np.squeeze(eps_data[:,:,round(cellSizeZ/2)])

ez_data = sim.get_array(center=mp.Vector3(),size=cell,component=mp.Ez)
#ez_data_mod = np.squeeze(ez_data[:,:,round(cellSizeZ/2)])
# ------------------------------------------------------------------ #
# Visualize geometry
# ------------------------------------------------------------------ #

plt.figure(dpi=150)
plt.imshow(eps_data.transpose(),interpolation='spline36',cmap='binary')
plt.imshow(np.real(ez_data.transpose()),interpolation='spline36',cmap='RdBu', alpha=0.9)
plt.axis('off')
plt.show()
