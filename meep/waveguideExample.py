# ---------------------------------------------------------------------------- #
# waveguideExample.py
# ---------------------------------------------------------------------------- #
#
#
# ---------------------------------------------------------------------------- #
# Import libraries
# ---------------------------------------------------------------------------- #

from __future__ import division

# Get project library path to import library files
import sys
import os
d = os.path.dirname(os.getcwd())
libPath = os.path.abspath(os.path.join(d, 'lib'))
sys.path.insert(0, libPath)

import meep as mp
import SiP_Materials
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------- #
# Setup Constants
# ---------------------------------------------------------------------------- #

# Speed of light in vacuum, m/s
c = 299792458

# default unit length is 1 um
um_scale = 1.0

# import silicon
Si   = mp.Medium(index=3.4757)#SiP_Materials.Si

# import SiO2
SiO2 = mp.Medium(index=1.4440)#SiP_Materials.SiO2

# ---------------------------------------------------------------------------- #
# Draw geometries
# ---------------------------------------------------------------------------- #

# Define slab waveguide dimensions (ignore 3rd dimension)
width     = 0.5
height    = 0.22
thickness = 1e20
waveguide = mp.Block(
    size=mp.Vector3(width,height,thickness),
    center=mp.Vector3(0,0),
    material=Si)
geometry = [waveguide]

# ---------------------------------------------------------------------------- #
# Run meep to extract omega-k relationship
# ---------------------------------------------------------------------------- #
# Specify size of bondary layers
boundaryCondWidth = 1.0
pml_layers = [mp.PML(boundaryCondWidth)]

# Specify cell sizes
cellX = 2*width + 2*boundaryCondWidth
cellY = 2*height + 2*boundaryCondWidth
cell = mp.Vector3(cellX,cellY,0)
resolution = 30

# Specify a broadband source
fcen = 1/1.55
df = 1
#sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ez, center=mp.Vector3())]
sources = [mp.Source(mp.ContinuousSource(fcen), component=mp.Ez, center=mp.Vector3(0,0))]
# Specify which kpts to choose
lambdaSource = 1.550
lambdaMax = 1.550001
lambdaMin = 1.549999
kmin = 1/lambdaMax
kmax = 1/lambdaMin
kSource = 1/lambdaSource

kpts = [mp.Vector3(kmin), mp.Vector3(kSource), mp.Vector3(kmax)]

sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    default_material=SiO2,
    resolution=resolution,
    force_complex_fields=True,
    symmetries = []
)

sim.init_fields()
sim.solve_cw()
#sim.output_efield()


'''
# Extract out omega from each k vector
all_freqs = sim.run_k_points(400, kpts)  # a list of lists of frequencies

# Get wavelength data in microns
k = np.asarray([v.x for v in kpts])

# Get Omega
f = np.asarray(all_freqs)


# get Wavelength
lambdaEst = 1 / f[:,0]

# Get effective index info
n_meep = k/f[:,0]


# Get group index info
ng = np.mean(np.diff(k) / np.diff(f))

print(f.size)
print(k)

'''


# ---------------------------------------------------------------------------- #
# Plot waveguide
# ---------------------------------------------------------------------------- #

eps_data = sim.get_array(center=mp.Vector3(),size=cell,component = mp.Dielectric)
ez_data = sim.get_array(center=mp.Vector3(),size=cell,component = mp.Ez)
plt.figure(1)
plt.imshow(np.abs(eps_data.transpose()),interpolation='spline36',cmap='binary')
plt.imshow(np.abs(ez_data.transpose()),interpolation='spline36',cmap='Reds',alpha=0.9)
plt.show()

plt.figure(2)
plt.imshow(np.abs(np.fft.fft2(ez_data.transpose())),interpolation='spline36',cmap='Reds',alpha=0.9)
plt.show()
# ---------------------------------------------------------------------------- #
# Plot refractive index data
# ---------------------------------------------------------------------------- #
'''
plt.figure(num=None, figsize=(6,5), dpi=100, facecolor='w', edgecolor='k')

plt.subplot(211)
plt.plot(lambdaPlot*1e6,np.sqrt(epsilonTheory))
plt.hold(True)
plt.plot(expdata[:,0],expdata[:,1],'o')
plt.plot(lambdaEst,n_meep,'o')
plt.hold(False)
plt.xlim(1,2)
plt.grid(True)

plt.ylabel('Index of Refraction')
plt.legend(('Analytic Fit','Paliks Handbook','Meep Simulation'))

plt.subplot(212)
plt.plot(lambdaEst,err,'o')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('$\%$ Error')
plt.grid(True)
plt.savefig('refractiveIndex.png')
plt.show()
'''
