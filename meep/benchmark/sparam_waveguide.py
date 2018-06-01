# ---------------------------------------------------------------------------- #
# sparam_waveguide.py
# ---------------------------------------------------------------------------- #
# Computes the S-Parameters of a simple Si waveguide and compares them to the
# corresponding analytical solution.
#
# TODO:
#
# ---------------------------------------------------------------------------- #
# VERSION
# ---------------------------------------------------------------------------- #
# 31 May 2018 - AMH - Initialization
#
# ---------------------------------------------------------------------------- #
# Import Libraries
# ---------------------------------------------------------------------------- #

import meep as mp
from meep import mpb
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import scipy.io as sio

# ---------------------------------------------------------------------------- #
# Important constants
# ---------------------------------------------------------------------------- #

# Debug parameters
debugSim = False
drawSimpleWaveguide = False

# Physical constants
eps0 = 8.854e-12;              # Permittivity of free space
mu0  = 4 * np.pi * 1e-7;       # Permeability of free space
c0   = 1/np.sqrt(mu0*eps0);    # Speed of light in free space (m/s)

# Material definitions
epsSi = 12
nSi   = np.sqrt(epsSi)
Si    = mp.Medium(epsilon=epsSi)

# Geometry definitions
waveguideWidth   = 0.5   # width of waveguide (y direction)
waveguideLength  = mp.inf    # length of each waveguide (x direction)
portLength       = 10    # gap between ports

# ---------------------------------------------------------------------------- #
# Simulate in Meep
# ---------------------------------------------------------------------------- #

# Create geometry
waveguide = mp.Block(
    material = Si,
    center=mp.Vector3(0,0,0),
    size=mp.Vector3(waveguideLength,waveguideWidth,mp.inf)
)

geometry = [waveguide]

# computational grid resolution
resolution = 40

# Formulate boundary layers
dpml = 2  # size of boundary layer
boundary_layers = [ mp.PML(dpml) ]

# Create computational cell
sx = 2*dpml + 2*portLength        # computational grid width
sy = 2*dpml + 2*waveguideWidth    # computational grid height
sz = 0 # computational grid thickness
cell_size = mp.Vector3(sx, sy, sz)  # computational grid

# Generate source
fLower = 1/2                    # lower frequency
fUpper = 1/1                    # upper frequency
fcen = np.mean([fLower,fUpper]) # center frequency
fwidth = (fUpper - fLower)        # frequency width
xsrc = -portLength/2
sources = [mp.EigenModeSource(mp.GaussianSource(frequency=fcen,width=fwidth),
                    direction = mp.X,
                    center=mp.Vector3(xsrc,0))]
# Simulation block
sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=boundary_layers,
                    geometry=geometry,
                    sources=sources)

# Number of frequency bins to compute
nfreq = 50

# Height of each monitor
monitorHeight = sy

# Port 1 Monitor
xm1 = -portLength/2 + 1  # x-coordinate of monitor
mflux1 = sim.add_eigenmode(fcen, fwidth, nfreq,
                          mp.FluxRegion(center=mp.Vector3(xm1,0,0),
                                        size=mp.Vector3(0,monitorHeight,mp.inf),
                                        direction=mp.X))

# Port 2 Monitor
xm2 = portLength/2 # x-coordinate of monitor
mflux2 = sim.add_eigenmode(fcen, fwidth, nfreq,
                          mp.FluxRegion(center=mp.Vector3(xm2,0,0),
                                        size=mp.Vector3(0,monitorHeight,mp.inf),
                                        direction=mp.X))
# Run the simulation
if debugSim:
    sim.run(until=5)
else:
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(0,0), 1e-8))

# Extract the S Parameters
bands = [1]  # indices of modes for which to compute expansion coefficients
f1 = sim.get_eigenmode_coefficients(mflux1, bands) # Fluxes on port 1
f2 = sim.get_eigenmode_coefficients(mflux2, bands) # Fluxes on port 2

a1 = f1[0,:,0] # Forward propogating wave on port 1
b1 = f1[0,:,1] # Backward propogating wave on port 1

a2 = f2[0,:,0] # Forward propogating wave on port 2
b2 = f2[0,:,1] # Backward propogating wave on port 2

S11 = b1 / a1
S12 = a2 / a1

# Pull corresponding frequency data for plots
freqSim = np.array(mp.get_flux_freqs(mflux1))
freqSimAdjust = freqSim * 1e6 * c0 * 1e-12

# Pull effective indices and k vector for each omega point
geometry_lattice = mp.Lattice(size=mp.Vector3(0, sy, 0))  # computational grid
ms = mpb.ModeSolver(
    geometry_lattice=geometry_lattice,
    geometry=geometry,
    resolution=resolution
)

findOmega = freqSim        # omegas for which we wish to find corresponding k's
k_meep = np.zeros(nfreq)   # preallocate k list
for iter in range(0,nfreq):
    values = ms.find_k(mp.NO_PARITY,                # polarization of interest
                       findOmega[iter],             # omega corresponding to k
                       1,                           # number of bands (min)
                       1,                           # number of bands (max)
                       mp.Vector3(1),               # direction in K space
                       1e-3,                        # convergence tolerance
                       findOmega[iter] * nSi,       # Guess for K
                       findOmega[iter] * 0.1,       # Min magnitude for K
                       findOmega[iter] * (nSi + 2)  # Max magnitude for K
                       )

    # Extract found k vector from list and store in vector
    k_meep[iter] = values[0]

# From found k vectors, calculate effective index and store
neff = k_meep / freqSim

# Save raw S Param to matlab file for debugging
sio.savemat('sparam_waveguide_data.mat', mdict={'freqs':freqSim,'a1': a1,'a2': a2,'b1': b1,'b2': b2})

# ---------------------------------------------------------------------------- #
# Generate analytical waveguide S parameters (magnitude and phase)
# ---------------------------------------------------------------------------- #

# Interpolate k vectors
numFreqPts = 1e3
freqAn     = np.linspace(fLower,fUpper,numFreqPts)
f          = interpolate.PchipInterpolator(freqSim, k_meep)
kInterp    = f(freqAn)

# Find distance between the ports
L = (xm2 - xm1)

# Generate the S12/S11 parameters
S12An = np.exp(1j*2*np.pi*kInterp*L)
S11An = np.zeros(int(numFreqPts))

# Convert frequency to THz
freqAn     = c0*freqAn *1e6*1e-12

# ---------------------------------------------------------------------------- #
# Plot and compare
# ---------------------------------------------------------------------------- #

plt.figure()

plt.subplot(311)
eps_data = sim.get_array(center=mp.Vector3(), size=cell_size, component=mp.Dielectric)
plt.imshow(np.rot90(eps_data), interpolation='spline36', cmap='binary',
           extent=[-sx/2,sx/2,-sy/2,sy/2])
plt.plot([xm1,xm1],[monitorHeight/2,-monitorHeight/2],color='r')
plt.plot([xm2,xm2],[monitorHeight/2,-monitorHeight/2],color='b')
plt.scatter([xsrc],[0],color='g')
plt.legend(('M1','M2','SRC'))
plt.xlabel('X ($\mu m$)')
plt.ylabel('Y ($\mu m$)')

plt.subplot(312)
plt.plot(freqAn,np.abs(S12An))
plt.plot(freqSimAdjust,np.abs(S12)**2,'.')
plt.plot(freqAn,np.abs(S11An))
plt.plot(freqSimAdjust,np.abs(S11)**2,'.')
plt.legend(('S12 - analytic','S12 - numeric','S11 - analytic','S11 - numeric'))
plt.xlabel('Frequency (THz)')
plt.ylabel('Magnitude')

plt.subplot(313)
angleAn = np.unwrap(np.angle(S12An))
plt.plot(freqAn,angleAn)
plt.plot(freqSimAdjust,np.unwrap(np.angle(S12)),'.')
plt.xlabel('Frequency (THz)')
plt.ylabel('phase (rad)')
plt.legend(('S12 - analytic','S12 - numeric'))

plt.tight_layout()

plt.show()
