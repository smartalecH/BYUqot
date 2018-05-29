# ---------------------------------------------------------------------------- #
# sparam.py
# ---------------------------------------------------------------------------- #
# A simple script that illustrates how to extract s parameters from a simple
# waveguide using meep.
#
# TODO:
#
# ---------------------------------------------------------------------------- #
# VERSION
# ---------------------------------------------------------------------------- #
# 29 May 2018 - AMH - Initialization
#
# ---------------------------------------------------------------------------- #
# Import Libraries
# ---------------------------------------------------------------------------- #

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------- #
# Important constants and materials
# ---------------------------------------------------------------------------- #

eps0 = 7.40523989e+02
susceptibilities = [mp.DrudeSusceptibility(frequency=5.28366673e-01, gamma=1.62316031, sigma=1.17034609e-03),
                    mp.LorentzianSusceptibility(frequency=4.19962269e+02, gamma=1.84154007e+03, sigma=7.44917949e+02)]
FR4 = mp.Medium(epsilon=eps0, E_susceptibilities=susceptibilities)

# ---------------------------------------------------------------------------- #
# Create Geometry
# ---------------------------------------------------------------------------- #
ww = 0.5      # waveguide width (microns)
wl = 10       # waveguie height (microns)

waveguide = mp.Block(
    material = FR4,
    center=mp.Vector3(0,0,0),
    size=mp.Vector3(wl,ww,mp.inf)
)

geometry = [waveguide]

# ---------------------------------------------------------------------------- #
# Setup Simulation
# ---------------------------------------------------------------------------- #

dpml = 3.0                         # size of boundary layer
sx = wl              # computational grid width
sy = dpml + 2*ww + dpml            # computaitonal grid height
cell_size = mp.Vector3(sx, sy, 0)  # computational grid

boundary_layers = [ mp.PML(dpml) ] # boundary layer

resolution = 10                    # computational grid resolution

lamcen = 1.55         # mode wavelength (microns)
fcen   = 1/lamcen     # mode frequency
sources = [ mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=0.5*fcen),
                               component=mp.Ez,
                               size=mp.Vector3(0,sy-2*dpml,0),
                               center=mp.Vector3(-0.5*sx+dpml,0,0),
                               eig_match_freq=True,
                               eig_parity=mp.ODD_Z,
                               eig_kpoint=mp.Vector3(0.4,0,0),
                               eig_resolution=32) ]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=boundary_layers,
                    geometry=geometry,
                    sources=sources)

# ---------------------------------------------------------------------------- #
# Extract S Parameters
# ---------------------------------------------------------------------------- #

xm = 0.5*sx - dpml  # x-coordinate of monitor
mflux = sim.add_eigenmode(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(xm,0), size=mp.Vector3(0,sy-2*dpml)))

sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(xm,0), 1e-10))

bands = [1]  # indices of modes for which to compute expansion coefficients
alpha = sim.get_eigenmode_coefficients(mflux, bands)

alpha0Plus  = alpha[2*0 + 0]  # coefficient of forward-traveling fundamental mode
alpha0Minus = alpha[2*0 + 1]  # coefficient of backward-traveling fundamental mode

print(alpha0Plus)
print(alpha0Minus)

# ---------------------------------------------------------------------------- #
# Extract Analytical S Parameters
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Plot and Compare
# ---------------------------------------------------------------------------- #
eps_data = sim.get_array(center=mp.Vector3(), size=cell_size, component=mp.Dielectric)
plt.figure(dpi=100)
plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
plt.axis('off')





plt.show()
