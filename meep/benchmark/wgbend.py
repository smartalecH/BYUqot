# ------------------------------------------------------------------------------ #
# From the Meep tutorial: transmission around a 90-degree waveguide bend in 2d.
# ------------------------------------------------------------------------------ #

from __future__ import division

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------ #
# Setup simulation domain
# ------------------------------------------------------------------------------ #

debugSim = 2

resolution = 10 # pixels/um

sx = 16  # size of cell in X direction
sy = 32  # size of cell in Y direction
cell = mp.Vector3(sx,sy,0)

dpml = 1.0
pml_layers = [mp.PML(dpml)]

pad = 4  # padding distance between waveguide and cell edge
w = 1    # width of waveguide

wvg_xcen =  0.5*(sx-w-2*pad)  # x center of vert. wvg
wvg_ycen = -0.5*(sy-w-2*pad)  # y center of horiz. wvg

monitorSize = 2*w

# ------------------------------------------------------------------------------ #
# Simulate straight waveguide using flux method
# ------------------------------------------------------------------------------ #

geometry = [mp.Block(size=mp.Vector3(mp.inf,w,mp.inf),
                     center=mp.Vector3(0,wvg_ycen,0),
                     material=mp.Medium(epsilon=12))]

fcen = 0.15  # pulse center frequency
df = 0.1     # pulse width (in frequency)

xsrc = -0.5*sx+dpml
ysrc = wvg_ycen
sources = [mp.EigenModeSource(mp.GaussianSource(fcen,fwidth=df), component=mp.Ez,
                     eig_match_freq=True, eig_band = 1,
                     center=mp.Vector3(xsrc,ysrc,0), size=mp.Vector3(0,2*monitorSize,0))]

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    resolution=resolution)

nfreq = 100  # number of frequencies at which to compute flux

# reflected flux
xm1 = xsrc+6
ym1 = wvg_ycen
refl_fr = mp.FluxRegion(center=mp.Vector3(xm1,ym1,0), size=mp.Vector3(0,monitorSize,0), direction=mp.X)
refl    = sim.add_eigenmode(fcen, df, nfreq, refl_fr)

# transmitted flux
xm2 = 0.5*sx-dpml
ym2 = wvg_ycen
tran_fr = mp.FluxRegion(center=mp.Vector3(xm2,ym2,0), size=mp.Vector3(0,monitorSize,0), direction=mp.X)
tran    = sim.add_eigenmode(fcen, df, nfreq, tran_fr)

pt = mp.Vector3(0.5*sx-dpml-0.5,wvg_ycen)

sim.run(until_after_sources=mp.stop_when_fields_decayed(50,mp.Ez,pt,1e-10))

# Pull the geometry of the first simulation
eps_data_wg = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)

# Eigenmode method results
bands = [1]
f1_simple = sim.get_eigenmode_coefficients(refl, bands) # Fluxes on port 1
f2_simple = sim.get_eigenmode_coefficients(tran, bands) # Fluxes on port 2
modeNumber = 0
a1_simple = f1_simple[modeNumber,:,0] # Forward propogating wave on port 1
b1_simple = f1_simple[modeNumber,:,1] # Backward propogating wave on port 1
a2_simple = f2_simple[modeNumber,:,0] # Forward propogating wave on port 1
b2_simple = f2_simple[modeNumber,:,1] # Backward propogating wave on port 1

flux_freqs     = np.array(mp.get_flux_freqs(refl))
wl = 1/flux_freqs
# ------------------------------------------------------------------------------ #
# Plot everything
# ------------------------------------------------------------------------------ #
plt.figure()
plt.plot(wl,abs(a1_simple)**2,label='a1')
plt.plot(wl,abs(b1_simple)**2,label='b1')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=4, mode="expand", borderaxespad=0)

plt.savefig('wgbend_1.png')

plt.figure()
plt.plot(wl,abs(a2_simple)**2,label='a2')
plt.plot(wl,abs(b2_simple)**2,label='b2')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=4, mode="expand", borderaxespad=0)

plt.savefig('wgbend_2.png')
