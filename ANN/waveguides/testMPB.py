
# ------------------------------------------------------------------------------ #
# Import libraries
# ------------------------------------------------------------------------------ #

import meep as mp
from   meep import mpb
import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys

# ------------------------------------------------------------------------------ #
# Setup simulation domain
# ------------------------------------------------------------------------------ #

c = 299792458    # Speed of light in vacuum, m/s

# Silicon dispersion relation
eps0Si = 7.98737492
epsLorentzSi = 3.68799143
omega0Si = 3.93282466e15
epsSi = lambda lam: eps0Si + epsLorentzSi*omega0Si ** 2 / (omega0Si ** 2 - (2*np.pi*c*1e6/lam)**2)

# Silicon Dioxide dispersion relation
eps0SiO2 = 2.119881
epsLorentzSiO2 = 49.43721
omega0SiO2 = 3.309238e13
epsSiO2 = lambda lam: eps0SiO2 + epsLorentzSiO2*omega0SiO2 ** 2 / (omega0SiO2 ** 2 - (2*np.pi*c*1e6/lam)**2)

sc_y = 2  # supercell width (um)
sc_z = 2  # supercell height (um)
geometry_lattice = mp.Lattice(size=mp.Vector3(0, sc_y, sc_z))
resolution = 32  # pixels/um

ms = mpb.ModeSolver(
    geometry_lattice=geometry_lattice,
    resolution=resolution,
    ensure_periodicity=False
)

# ------------------------------------------------------------------------------ #
# Loop through experiment
# ------------------------------------------------------------------------------ #

# get current frequency and wavelength
currentLambda = 1.55
omegaIn       = 1 / currentLambda

# Update material parameters
nSi = epsSi(currentLambda)
Si = mp.Medium(epsilon=nSi)

nSiO2 = epsSiO2(currentLambda)
SiO2 = mp.Medium(epsilon=nSiO2)

# Update geometry
wCurrent = 0.5
hCurrent = 0.22
ms.geometry = [mp.Block(size=mp.Vector3(mp.inf, wCurrent, hCurrent), material=Si)]
ms.default_material = SiO2

# Find TE modes vector
kTE = np.array(
    ms.find_k(
        p              = mp.EVEN_Z,       # Polarization
        omega          = omegaIn,            # Omega to find corresponding k
        band_min       = 1,                  # Minimum band index
        band_max       = 2,                  # Max band index
        korig_and_kdir = mp.Vector3(1),      # K vector orientation
        tol            = 1e-3,               # solver tolerance
        kmag_guess     = omegaIn * 3.45,     # initial guess
        kmag_min       = omegaIn * 0.1,      # Minimum
        kmag_max       = omegaIn * 4))       # Maximum

# Find TM modes
kTM = np.array(
    ms.find_k(
        p              = mp.ODD_Z,       # Polarization
        omega          = omegaIn,            # Omega to find corresponding k
        band_min       = 1,                  # Minimum band index
        band_max       = 2,                  # Max band index
        korig_and_kdir = mp.Vector3(1),      # K vector orientation
        tol            = 1e-3,               # solver tolerance
        kmag_guess     = omegaIn * 3.45,     # initial guess
        kmag_min       = omegaIn * 0.1,      # Minimum
        kmag_max       = omegaIn * 4))       # Maximum

print(kTE/omegaIn)
print(kTM/omegaIn)
