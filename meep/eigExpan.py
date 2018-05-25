# ---------------------------------------------------------------------------- #
# Import libraries
# ---------------------------------------------------------------------------- #

import meep as mp
from meep import mpb
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------- #
# Materials, constants, and common cimulation params
# ---------------------------------------------------------------------------- #

Si = mp.Medium(index=3.4757)
SiO2 = mp.Medium(index=1.444)

ww      = 0.5   # waveguide width (microns)
wh      = 0.22  # waveguide height (microns)
rradius = 5     # ring radius (microns)
gap     = .1   # gap between ring and bus (microns)
conLen  = 5     # connector length (microns)

resolution = 10

PML = 2.0
pml_layers = [mp.PML(PML)]

cellWidth = 3*rradius + 2*PML
cellHeight = 3*rradius + 2*PML

cell = mp.Vector3(cellWidth,cellHeight,0)

force_complex_fields = True

fcen = 1/1.55

xm = 1.25*rradius  # x-coordinate of monitor
ym = rradius + ww + gap # y-coordinate of monitor

# ---------------------------------------------------------------------------- #
# Simulate Simple Waveguide
# ---------------------------------------------------------------------------- #

# ------------ Geometries ------------#
bus =  mp.Block(
    size=mp.Vector3(1e20,ww,wh),
    center = mp.Vector3(0,rradius + gap + ww,0),
    material = Si
)
geometry = [bus]

# ------------ Simulation ------------#

sources = [ mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=0.5*fcen),
                               #component=mp.Ez,
                               size=mp.Vector3(0,ww,0),
                               center=mp.Vector3(-1.25*rradius,rradius + ww + gap,0),
                               #center = mp.Vector3(rradius,-rradius,0),
                               eig_match_freq=True,
                               eig_parity=mp.ODD_Z,
                               eig_kpoint=mp.Vector3(0.4,0,0),
                               eig_resolution=32,
                               direction = mp.X
                               )]


sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=pml_layers,
    force_complex_fields=force_complex_fields,
    resolution=resolution,
    default_material = SiO2
)


fluxOrig = sim.add_eigenmode(fcen, 0.5*fcen, 100, mp.FluxRegion(center=mp.Vector3(xm,ym), size=mp.Vector3(0,2*ww)))
sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(xm,ym), 1e-10))

rawFreqs = mp.get_flux_freqs(fluxOrig)
rawFluxes = mp.get_fluxes(fluxOrig)

# ---------------------------------------------------------------------------- #
# Simulate half ring
# ---------------------------------------------------------------------------- #

outerRing = mp.Cylinder(
    center = mp.Vector3(0,0,0),
    radius = rradius + ww/2,
    height = wh,
    material = Si
)

innerRing = mp.Cylinder(
    center = mp.Vector3(0,0,0),
    radius = rradius - ww/2,
    height = wh,
    material = SiO2
)

bus =  mp.Block(
    size=mp.Vector3(1e20,ww,wh),
    center = mp.Vector3(0,rradius + gap + ww,0),
    material = Si
)

ringBlock = mp.Block(
    size = mp.Vector3(1e20,4*rradius,1e20),
    center = mp.Vector3(0,-2*rradius,0),
    material = SiO2
)

ringConnector1 = mp.Block(
    size = mp.Vector3(ww,6*rradius,wh),
    center = mp.Vector3(rradius,-3*rradius,0),
    material = Si
)
ringConnector2 = mp.Block(
    size = mp.Vector3(ww,6*rradius,wh),
    center = mp.Vector3(-rradius,-3*rradius,0),
    material = Si
)

geometry = [outerRing,innerRing,bus,ringBlock,ringConnector1,ringConnector2]
#geometry = [outerRing,innerRing,bus]

sources = [ mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=0.5*fcen),
                               #component=mp.Ez,
                               size=mp.Vector3(0,ww,0),
                               center=mp.Vector3(-1.25*rradius,rradius + ww + gap,0),
                               #center = mp.Vector3(rradius,-rradius,0),
                               eig_match_freq=True,
                               eig_parity=mp.ODD_Z,
                               eig_kpoint=mp.Vector3(0.4,0,0),
                               eig_resolution=32,
                               direction = mp.X
                               )]

sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=pml_layers,
    force_complex_fields=force_complex_fields,
    resolution=resolution,
    default_material = SiO2
)

xm = 1.25*rradius  # x-coordinate of monitor
ym = rradius + ww + gap # y-coordinate of monitor
fluxOrig = sim.add_eigenmode(fcen, 0.5*fcen, 100, mp.FluxRegion(center=mp.Vector3(xm,ym), size=mp.Vector3(0,2*ww)))
sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(xm,ym), 1e-10))

postFreqs = mp.get_flux_freqs(fluxOrig)
postFluxes = mp.get_fluxes(fluxOrig)

response = np.array(postFluxes) / np.array(rawFluxes)
lam = 1/np.array(postFreqs)
# ---------------------------------------------------------------------------- #
# Plot everything
# ---------------------------------------------------------------------------- #

eps_data = sim.get_array(center=mp.Vector3(),size=cell,component = mp.Dielectric)
ez_data = sim.get_array(center=mp.Vector3(),size=cell,component = mp.Ez)

plt.subplot(211)
plt.imshow(np.abs(eps_data.transpose()),interpolation='spline36',cmap='binary')

plt.subplot(212)
plt.plot(lam,response)
plt.show()
