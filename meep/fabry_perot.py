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
from meep import mpb
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

# ---------------------------------------------------------------------------- #
# Important constants and materials
# ---------------------------------------------------------------------------- #

# Debug parameters
debugSim = False
drawSimpleWaveguide = False

# Physical constants
eps0 = 8.854e-12;           # Permittivity of free space
mu0  = 4 * np.pi * 1e-7;       # Permeability of free space
c0   = 1/np.sqrt(mu0*eps0);    # Speed of light in free space (m/s)

# Material definitions
epsSi = 12
epsGlass = 2.07
nSi   = np.sqrt(epsSi)
n0    = 1
Air   = mp.Medium(epsilon=1)
Si    = mp.Medium(epsilon=epsSi)
SiO2  = mp.Medium(epsilon=epsGlass)

ngap = np.sqrt(1)

# Geometry definitions
waveguideWidth = 0.5    # width of waveguide (y direction)
waveguideLength = 30    # length of each waveguide (x direction)
gap             = 10    # gap between waveguides


################################################################################
################################################################################
# Draw simple waveguide for simple calculations
################################################################################
################################################################################
if drawSimpleWaveguide:

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
    mflux1 = sim.add_flux(fcen, 2*fwidth, nfreq,
                              mp.FluxRegion(center=mp.Vector3(xm1,0,0),
                                            size=mp.Vector3(0,monitorHeight,mp.inf),
                                            direction=mp.X))

    # ------------------- Port 2 Monitor -------------------- #

    xm2 = sx/4 # x-coordinate of monitor
    mflux2 = sim.add_flux(fcen, 2*fwidth, nfreq,
                              mp.FluxRegion(center=mp.Vector3(xm2,0,0),
                                            size=mp.Vector3(0,monitorHeight,mp.inf),
                                            direction=mp.X))

    # ------------------- Run Simulation -------------------- #

    if debugSim:
        sim.run(until=5)
    else:
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
# Extract effective index of cavity and view analytical results
################################################################################
################################################################################

gapGeo = mp.Block(
    material = SiO2,
    center=mp.Vector3(0,0,0),
    size=mp.Vector3(mp.inf,waveguideWidth,mp.inf)
)

geometry = [gapGeo]

sy = 2 * waveguideWidth  # computational grid height
geometry_lattice = mp.Lattice(size=mp.Vector3(0, sy, 0))  # computational grid

# The k (i.e. beta, i.e. propagation constant) points to look at, in
# units of 2*pi/um.  We'll look at num_k points from k_min to k_max.
num_k = 9
k_min = 0.1
k_max = 3.0
k_points = [mp.Vector3(k_min), mp.Vector3(k_max)]

resolution = 32  # pixels/um

# Increase this to see more modes.  (The guided ones are the ones below the
# light line, i.e. those with frequencies < kmag / 1.45, where kmag
# is the corresponding column in the output if you grep for "freqs:".)
num_bands = 1

ms = mpb.ModeSolver(
    geometry_lattice=geometry_lattice,
    geometry=geometry,
    k_points=k_points,
    resolution=resolution,
    num_bands=num_bands
)

numWavelength = 250
wavelengthVec = np.linspace(1,2,numWavelength)
findOmega = 1 / wavelengthVec  # frequency corresponding to wavelengths of interest

# Output the x component of the Poynting vector for num_bands bands at omega
k_meep = np.zeros(numWavelength)
for iter in range(0,numWavelength):
    values = ms.find_k(mp.NO_PARITY, findOmega[iter], 1, 1, mp.Vector3(1), 1e-3, findOmega[iter] * 3.45,
              findOmega[iter] * 0.1, findOmega[iter] * 4)
    k_meep[iter] = values[0]


neff = k_meep * wavelengthVec

# Solve for the analytical structure
L = gap
R  = ((nSi - n0) / (nSi + n0)) ** 2
FP = lambda k: (1-R) ** 2 / ((1-R) ** 2 + 4 * R * (np.sin(k*L)) ** 2)
k = 2*np.pi*k_meep
T = FP(k)
R = 1-T
fadjust = c0 / (wavelengthVec * 1e-6) * 1e-12

################################################################################
################################################################################
# Draw fabry perot waveguide for actual calculations
################################################################################
################################################################################

# ------------------- Draw Geometry -------------------- #

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

resolution = 40 # computational grid resolution

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

# ------------------- Post-processing -------------------- #

nfreq = 500 # Number of frequency bins to compute
monitorHeight = sy

# ------------------- Port 1 Monitor -------------------- #

xm1 = -sx/4 + 1/resolution  # x-coordinate of monitor
mflux1 = sim.add_eigenmode(fcen, 2*fwidth, nfreq,
                          mp.FluxRegion(center=mp.Vector3(xm1,0,0),
                                        size=mp.Vector3(0,monitorHeight,mp.inf),
                                        direction=mp.X))

# ------------------- Port 2 Monitor -------------------- #

xm2 = sx/4 # x-coordinate of monitor
mflux2 = sim.add_eigenmode(fcen, 2*fwidth, nfreq,
                          mp.FluxRegion(center=mp.Vector3(xm2,0,0),
                                        size=mp.Vector3(0,monitorHeight,mp.inf),
                                        direction=mp.X))

# ------------------- Run Simulation -------------------- #

if debugSim:
    sim.run(until=5)
else:
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(0,0), 1e-8))

# ------------------- Extract S Params -------------------- #

'''
f1 = np.array(mp.get_fluxes(mflux1))
f2 = np.array(mp.get_fluxes(mflux2))

S12 = f1 / f1_simple
#sio.savemat('fabry_perot.mat', mdict={'freqs':freqs,'f1': f1,'f2': f2})
'''

bands = [1]  # indices of modes for which to compute expansion coefficients
f1 = sim.get_eigenmode_coefficients(mflux1, bands)
f2 = sim.get_eigenmode_coefficients(mflux2, bands)

a1 = f1[0,:,0]
b1 = f1[0,:,1]

a2 = f2[0,:,0]
b2 = f2[0,:,1]

S11 = b1 / a1
S12 = a2 / a1



freqs = np.array(mp.get_flux_freqs(mflux1))
waveSim = 1 / (freqs)

freqsSimAdjust = freqs * 1e6 * c0 * 1e-12

sio.savemat('fabry_perot.mat', mdict={'freqs':freqs,'a1': a1,'a2': a2,'b1': b1,'b2': b2})

# ------------------- Plot the waveguide structure -------------------- #
plt.figure(2)
plt.subplot(311)

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
