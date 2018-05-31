# ---------------------------------------------------------------------------- #
# fabry_perot.py
# ---------------------------------------------------------------------------- #
# A simple script that simulates a 2D Fabry-Perot cavity, extracts the
# S parameters, and compares them to the analytical formula.
#
# TODO:
#
# ---------------------------------------------------------------------------- #
# VERSION
# ---------------------------------------------------------------------------- #
# 30 May 2018 - AMH - Initialization
#
# ---------------------------------------------------------------------------- #
# Import Libraries
# ---------------------------------------------------------------------------- #

import meep as mp
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

# ---------------------------------------------------------------------------- #
# Important constants and materials
# ---------------------------------------------------------------------------- #
# Physical constants
eps0 = 8.854e-12;           # Permittivity of free space
mu0  = 4 * np.pi * 1e-7;       # Permeability of free space
c0   = 1/np.sqrt(mu0*eps0);    # Speed of light in free space (m/s)

epsSi = 12
epsGlass = 2.07
nSi   = np.sqrt(epsSi)
n0    = 1
Air   = mp.Medium(epsilon=1)
Si    = mp.Medium(epsilon=epsSi)
SiO2  = mp.Medium(epsilon=epsGlass)

ngap = np.sqrt(epsGlass)

################################################################################
################################################################################
# Draw simple waveguide for simple calculations
################################################################################
################################################################################

waveguideWidth = 0.5    # width of waveguide (y direction)
waveguideLength = 30    # length of each waveguide (x direction)
gap             = 10    # gap between waveguides

waveguide = mp.Block(
    material = Si,
    center=mp.Vector3(0,0,0),
    size=mp.Vector3(2*waveguideLength + gap,waveguideWidth,mp.inf)
)

geometry = [waveguide]

resolution = 20 # computational grid resolution

# ------------------- PML boundary layers -------------------- #

dpml = 1  # size of boundary layer
boundary_layers = [ mp.PML(dpml) ]

# ------------------- Computational cell size -------------------- #

sx = waveguideLength   # computational grid width
sy = 2*dpml + 2 * waveguideWidth  # computational grid height
sz = 0 # computational grid thickness
cell_size = mp.Vector3(sx, sy, sz)  # computational grid

# ------------------- Source -------------------- #

fLower = 1/2                    # lower frequency
fUpper = 1/1                    # upper frequency
fcen = np.mean([fLower,fUpper]) # center frequency
fwidth = (fUpper - fLower)        # frequency width

xsrc = -sx/4
sources = [mp.EigenModeSource(mp.GaussianSource(frequency=fcen,width=fwidth),
                    direction = mp.X,
                    center=mp.Vector3(xsrc,0))]

# ------------------- Simulation Block -------------------- #

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=boundary_layers,
                    geometry=geometry,
                    sources=sources)

nfreq = 200 # Number of frequency bins to compute
monitorHeight = sy-2*dpml

# ------------------- Port 1 Monitor -------------------- #

xm1 = -sx/4  # x-coordinate of monitor
mflux1 = sim.add_flux(fcen, fwidth, nfreq,
                          mp.FluxRegion(center=mp.Vector3(xm1,0,0),
                                        size=mp.Vector3(0,monitorHeight,mp.inf),
                                        direction=mp.X))

# ------------------- Port 2 Monitor -------------------- #

xm2 = sx/4 # x-coordinate of monitor
mflux2 = sim.add_flux(fcen, fwidth, nfreq,
                          mp.FluxRegion(center=mp.Vector3(xm2,0,0),
                                        size=mp.Vector3(0,monitorHeight,mp.inf),
                                        direction=mp.X))

# ------------------- Run Simulation -------------------- #

#sim.run(until=5)
sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(0,0), 1e-8))

# ------------------- Extract S Params -------------------- #
f1_simple = np.array(mp.get_fluxes(mflux1))
f2_simple = np.array(mp.get_fluxes(mflux2))

freqs = np.array(mp.get_flux_freqs(mflux1))
waveSim = 1 / (freqs)

freqsSimAdjust = freqs * 1e6 * c0 * 1e-12

eps_data = sim.get_array(center=mp.Vector3(), size=cell_size, component=mp.Dielectric)

plt.figure(1)
plt.subplot(311)
plt.imshow(np.rot90(eps_data), interpolation='spline36', cmap='binary',
           extent=[-sx/2,sx/2,-sy/2,sy/2])
plt.plot([xm1,xm1],[monitorHeight/2,-monitorHeight/2],color='r')
plt.plot([xm2,xm2],[monitorHeight/2,-monitorHeight/2],color='b')
plt.scatter([xsrc],[0],color='g')
plt.legend(('M1','M2','SRC'))
plt.xlabel('X ($\mu m$)')
plt.ylabel('Y ($\mu m$)')

plt.subplot(312)
plt.title('F1')
plt.plot(freqsSimAdjust,f1_simple)
plt.xlabel('Frequency (THz)')
plt.ylabel('Transmission')

plt.subplot(313)
plt.title('F2')
plt.plot(freqsSimAdjust,f1_simple)
#plt.plot(freqsSimAdjust,np.abs(S11) ** 2)
plt.xlabel('Frequency (THz)')
plt.ylabel('Transmission')

sim.reset_meep()

################################################################################
################################################################################
# Draw fabry perot waveguide for actual calculations
################################################################################
################################################################################

# ---------------------------------------------------------------------------- #
# Draw Geometry
# ---------------------------------------------------------------------------- #

###############################
# All units in terms of microns!
###############################

waveguideWidth = 0.5    # width of waveguide (y direction)
waveguideLength = 30    # length of each waveguide (x direction)
gap             = 10    # gap between waveguides

waveguide = mp.Block(
    material = Si,
    center=mp.Vector3(0,0,0),
    size=mp.Vector3(2*waveguideLength + gap,waveguideWidth,mp.inf)
)

gapGeo = mp.Block(
    material = SiO2,
    center=mp.Vector3(0,0,0),
    size=mp.Vector3(gap,waveguideWidth,mp.inf)
)

geometry = [waveguide,gapGeo]

# ---------------------------------------------------------------------------- #
# Setup Simulation
# ---------------------------------------------------------------------------- #

resolution = 20 # computational grid resolution

# ------------------- PML boundary layers -------------------- #

dpml = 2  # size of boundary layer
boundary_layers = [ mp.PML(dpml) ]

# ------------------- Computational cell size -------------------- #

sx = waveguideLength   # computational grid width
sy = 2*dpml + 2 * waveguideWidth  # computational grid height
sz = 0 # computational grid thickness
cell_size = mp.Vector3(sx, sy, sz)  # computational grid



# ------------------- Source -------------------- #

fLower = 1/2                    # lower frequency
fUpper = 1/1                    # upper frequency
fcen = np.mean([fLower,fUpper]) # center frequency
fwidth = (fUpper - fLower)        # frequency width

xsrc = -sx/4
sources = [mp.EigenModeSource(mp.GaussianSource(frequency=fcen,width=fwidth),
                    direction = mp.X,
                    center=mp.Vector3(xsrc,0))]

# ------------------- Simulation Block -------------------- #

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=boundary_layers,
                    geometry=geometry,
                    sources=sources)

# ---------------------------------------------------------------------------- #
# Post-processing
# ---------------------------------------------------------------------------- #

nfreq = 200 # Number of frequency bins to compute
monitorHeight = sy-2*dpml

# ------------------- Port 1 Monitor -------------------- #

xm1 = -sx/4  # x-coordinate of monitor
mflux1 = sim.add_flux(fcen, fwidth, nfreq,
                          mp.FluxRegion(center=mp.Vector3(xm1,0,0),
                                        size=mp.Vector3(0,monitorHeight,mp.inf),
                                        direction=mp.X))

# ------------------- Port 2 Monitor -------------------- #

xm2 = sx/4 # x-coordinate of monitor
mflux2 = sim.add_flux(fcen, fwidth, nfreq,
                          mp.FluxRegion(center=mp.Vector3(xm2,0,0),
                                        size=mp.Vector3(0,monitorHeight,mp.inf),
                                        direction=mp.X))

# ------------------- Run Simulation -------------------- #

#sim.run(until=5)
sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(0,0), 1e-8))

# ------------------- Extract S Params -------------------- #
f1 = np.array(mp.get_fluxes(mflux1))
f2 = np.array(mp.get_fluxes(mflux2))



S12 = f1 / f1_simple
'''
bands = [1]  # indices of modes for which to compute expansion coefficients
f1 = sim.get_eigenmode_coefficients(mflux1, bands)
f2 = sim.get_eigenmode_coefficients(mflux2, bands)
'''
'''
a1 = f1[0,:]
a2 = f2[0,:]
b1 = f1[1,:]
b2 = f2[1,:]

S11 = b1 / a1
S12 = a2 / a1
'''


freqs = np.array(mp.get_flux_freqs(mflux1))
waveSim = 1 / (freqs)

freqsSimAdjust = freqs * 1e6 * c0 * 1e-12

#sio.savemat('fabry_perot.mat', mdict={'freqs':freqs,'a1': a1,'a2': a2,'b1': b1,'b2': b2})
sio.savemat('fabry_perot.mat', mdict={'freqs':freqs,'f1': f1,'f2': f2})
# ---------------------------------------------------------------------------- #
# Analytical comparison
# ---------------------------------------------------------------------------- #

# Solve for the analytical structure
L = gap
R  = ((nSi - n0) / (nSi + n0)) ** 2

FP = lambda k: (1-R) ** 2 / ((1-R) ** 2 + 4 * R * (np.sin(k*L)) ** 2)
numWavelength = 1000
wavelengthVec = np.linspace(1,2,numWavelength)
k = 2*np.pi*ngap/wavelengthVec
T = FP(k)
R = 1-T

fadjust = c0 / (wavelengthVec * 1e-6) * 1e-12

plt.figure(2)
plt.subplot(311)

# ------------------- Plot the waveguide structure -------------------- #

eps_data = sim.get_array(center=mp.Vector3(), size=cell_size, component=mp.Dielectric)
ez_data = sim.get_array(center=mp.Vector3(), size=cell_size, component=mp.Ez)

plt.imshow(np.rot90(eps_data), interpolation='spline36', cmap='binary',
           extent=[-sx/2,sx/2,-sy/2,sy/2])
plt.plot([xm1,xm1],[monitorHeight/2,-monitorHeight/2],color='r')
plt.plot([xm2,xm2],[monitorHeight/2,-monitorHeight/2],color='b')
plt.scatter([xsrc],[0],color='g')
plt.legend(('M1','M2','SRC'))
plt.xlabel('X ($\mu m$)')
plt.ylabel('Y ($\mu m$)')

# ------------------- Plot Transmission -------------------- #

plt.subplot(312)
plt.title('Transmission ($|S_{12}|^2$)')
plt.plot(fadjust,T)
plt.plot(freqsSimAdjust,np.abs(S12) ** 2)
plt.xlabel('Frequency (THz)')
plt.ylabel('Transmission')


# ------------------- Plot Reflection -------------------- #
plt.subplot(313)
plt.title('Reflection ($|S_{11}|^2$)')
plt.plot(fadjust,R)
#plt.plot(freqsSimAdjust,np.abs(S11) ** 2)
plt.xlabel('Frequency (THz)')
plt.ylabel('Transmission')


plt.tight_layout()

plt.show()
