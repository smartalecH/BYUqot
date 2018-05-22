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
c = 3e8

# default unit length is 1 um
um_scale = 1.0

SiO2_range = mp.FreqRange(min=1/1.77, max=1/0.25)

SiO2_frq1 = 1/(0.103320160833333*um_scale)
SiO2_gam1 = 0#1/(12.3984193000000*um_scale)
SiO2_sig1 = 1.12

SiO2_susc = [ mp.LorentzianSusceptibility(frequency=SiO2_frq1, gamma=SiO2_gam1, sigma=SiO2_sig1) ]

SiO2 = mp.Medium(epsilon=1.0, E_susceptibilities=SiO2_susc, valid_freq_range=SiO2_range)

# ---------------------------------------------------------------------------- #
# Let's view silicon's theoretical refractive index plot
# ---------------------------------------------------------------------------- #

# Data from Lukas's handbook
eps = 1
eps_lorentz = SiO2_sig1
omega0 = 2*np.pi*3e8*SiO2_frq1*1e6
print(omega0)
delta0 = 0

# Generate data to plot over valid range
lambdaMin = 0.25e-6
lambdaMax = 1.77e-6
lambdaN   = 100
lambdaPlot = np.linspace(lambdaMin,lambdaMax,lambdaN)
num = eps_lorentz * (omega0**2)
den = (omega0**2 - 2*1j*delta0*2*np.pi * c / lambdaPlot - (2*np.pi*c/lambdaPlot)**2)
#epsilonTheory = eps +  num / den
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
# Let's simulate the w-k plot and eps plot
# ---------------------------------------------------------------------------- #

# Now let's rearrange the data to the form Meep uses
cell = mp.Vector3()
resolution = 30

fcen = 1
df = 4

sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ez, center=mp.Vector3())]

kmin = 1/1.7
kmax = 1/0.25
k_interp = 10

kpts = mp.interpolate(k_interp, [mp.Vector3(kmin), mp.Vector3(kmax)])

sim = mp.Simulation(
    cell_size=cell,
    geometry=[],
    sources=sources,
    default_material=SiO2,
    resolution=resolution
)

# Extract out omega from each k vector
all_freqs = sim.run_k_points(300, kpts)  # a list of lists of frequencies
print(all_freqs)

# Get wavelength data in microns
k = np.asarray([v.x for v in kpts])


# Get Omega
f         = np.asarray([v[0] for v in all_freqs])
lambdaEst = 1 / k


# Get refractive index info
n_meep = k/f

err = abs(np.sqrt(epsLam(lambdaEst)) - n_meep) / np.sqrt(epsLam(lambdaEst)) * 100
print('Error:')
print(err)



# ---------------------------------------------------------------------------- #
# Plot data
# ---------------------------------------------------------------------------- #

plt.plot(lambdaPlot*1e6,np.sqrt(epsilonTheory))
plt.hold(True)
#plt.plot(expdata[:,0],expdata[:,1],'o')
plt.plot(lambdaEst,n_meep,'o')
plt.hold(True)
plt.xlim(lambdaMin*1e6,lambdaMax*1e6)
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Index of Refraction')


plt.show()
