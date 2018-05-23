# Material dispersion example, from the Meep tutorial.  Here, we simply
# simulate homogenous space filled with a dispersive material, and compute
# its modes as a function of wavevector k.  Since omega/c = k/n, we can
# extract the dielectric function epsilon(omega) = (ck/omega)^2.
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

# Speed of light in vacuum, m/s
c = 299792458

# default unit length is 1 um
um_scale = 1.0

# ---------------------------------------------------------------------------- #
# Let's view silicon's theoretical refractive index plot
# ---------------------------------------------------------------------------- #

# Data from Lukas's handbook
eps = 7.9874
eps_lorentz = 3.6880
omega0 = 3.9328e15
delta0 = 0

# Generate data to plot over valid range
lambdaMin = 1e-6
lambdaMax = 2e-6
lambdaN   = 100
lambdaPlot = np.linspace(lambdaMin,lambdaMax,lambdaN)
num = eps_lorentz * (omega0**2)
den = (omega0**2 - 2*1j*delta0*2*np.pi * c / lambdaPlot - (2*np.pi*c/lambdaPlot)**2)
epsLam = lambda x: eps + (eps_lorentz * (omega0**2)) / ((omega0**2 - 2*1j*delta0*2*np.pi * c / (x*1e-6) - (2*np.pi*c/(x*1e-6))**2))
epsilonTheory = epsLam(lambdaPlot*1e6)

# ---------------------------------------------------------------------------- #
# Let's load experimental data
# ---------------------------------------------------------------------------- #

# Reference:
# H. H. Li. Refractive index of silicon and germanium and its wavelength and
# temperature derivatives, J. Phys. Chem. Ref. Data 9, 561-658 (1993)

filename = 'SiData.csv'
expdata = np.genfromtxt(fname=filename,delimiter=',',skip_header=1)

# ---------------------------------------------------------------------------- #
# Run meep to extract omega-k relationship
# ---------------------------------------------------------------------------- #

cell = mp.Vector3()
resolution = 200

fcen = 1
df = 3
sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ez, center=mp.Vector3())]

kmin = 1/0.8
kmax = 1/0.3
k_interp = 10

kpts = mp.interpolate(k_interp, [mp.Vector3(kmin), mp.Vector3(kmax)])

sim = mp.Simulation(
    cell_size=cell,
    geometry=[],
    sources=sources,
    default_material=SiP_Materials.Si,
    resolution=resolution,
)

# Extract out omega from each k vector
all_freqs = sim.run_k_points(300, kpts)  # a list of lists of frequencies

# Get wavelength data in microns
k = np.asarray([v.x for v in kpts])

# Get Omega
f = np.asarray([v[0].real for v in all_freqs])

# get Wavelength
lambdaEst = 1 / f

# Get refractive index info
n_meep = k/f

# Extract Error
y = np.sqrt(np.real(epsLam(lambdaEst)))
yhat = n_meep
err = 100 * abs((y - yhat)/y)

# ---------------------------------------------------------------------------- #
# Plot refractive index data
# ---------------------------------------------------------------------------- #

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
