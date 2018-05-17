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
ringRadius        = 10      # radius of the center of the ring waveguide
ringWaveguideWith = 0.5     # width of the ring waveguide
gapLength         = 0.2     # distance between ring edge and bus
topBusWidth       = 0.5     # width of top bus
bottomBusWidth    = 0.5     # width of bottom bus
topBus            = True    # Boolean - true if including top bus
bottomBus         = True    # Boolean - true if including bottom bus
racetrackLength   = 0       # length of the racetrack

waveguideThickness = 0.22   # thickness of all waveguides

# ------------------------------------------------------------------ #
# Draw geometry
# ------------------------------------------------------------------ #

# preallocate geometry
geometry = [];

# draw the ring
outerRingGeom = mp.Cylinder(
                    radius = ringRadius + ringWaveguideWith / 2,
                    height = waveguideThickness,
                    center = mp.Vector3(0,0,0),
                    material = Si
                    )
innerRingGeom = mp.Cylinder(
                    radius = ringRadius - ringWaveguideWith / 2,
                    height = waveguideThickness,
                    center = mp.Vector3(0,0,0),
                    material = SiO2
                    )
geometry.append(outerRingGeom)
geometry.append(innerRingGeom)

# draw the top bus
if topBus:
    yPos = ringRadius + ringWaveguideWith/2 + gapLength + topBusWidth/2;
    topBusGeom = mp.Block(
       mp.Vector3(2*ringRadius + ringWaveguideWith,bottomBusWidth,waveguideThickness),
       center = mp.Vector3(0,yPos,0),
       material = Si
       )
    geometry.append(topBusGeom)


# draw the bottom bus
if bottomBus:
    yPos = ringRadius + ringWaveguideWith/2 + gapLength + bottomBusWidth/2;
    bottomBusGeom = mp.Block(
       mp.Vector3(2*ringRadius + ringWaveguideWith,bottomBusWidth,waveguideThickness),
       center = mp.Vector3(0,-yPos,0),
       material = Si
       )
    geometry.append(bottomBusGeom)

# ------------------------------------------------------------------ #
# Setup simulation domain
# ------------------------------------------------------------------ #

# Size of domain
cell = mp.Vector3(2*ringRadius + ringWaveguideWith,2*ringRadius + ringWaveguideWith,3*waveguideThickness)

# Sources
sources = [mp.Source(mp.ContinuousSource(frequency=0.15),
            component=mp.Ez,
            center=mp.Vector3(0,0))]

# Boundary conditions
pml_layers = [mp.PML(1.0)]

# Resolution
resolution = 10

sim = mp.Simulation(cell_size = cell,
                    boundary_layers = pml_layers,
                    geometry = geometry,
                    sources = sources,
                    resolution = resolution)

eps_data = sim.get_array(center=mp.Vector3(),size=cell,component=mp.Dielectric)


# ------------------------------------------------------------------ #
# Visualize geometry
# ------------------------------------------------------------------ #
