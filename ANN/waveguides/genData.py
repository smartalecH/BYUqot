
# ------------------------------------------------------------------------------ #
# Import libraries
# ------------------------------------------------------------------------------ #

import meep as mp
from   meep import mpb
import matplotlib.pyplot as plt
import numpy as np
import h5py

# ------------------------------------------------------------------------------ #
# Define size of simulation array
# ------------------------------------------------------------------------------ #

nW      = int(2)
nH      = int(2)
nLambda = int(2)
lMin    = int(100)

# ------------------------------------------------------------------------------ #
# Setup simulation domain
# ------------------------------------------------------------------------------ #

# Silicon dispersion relations
c = 299792458    # Speed of light in vacuum, m/s
eps0 = 7.9874
epsLorentz = 3.6880
omega0 = 3.9328e15
epsSi = lambda lam: eps0 + epsLorentz*omega0 ** 2 / (omega0 ** 2 - (2*np.pi*c*1e6/lam)**2)

# Width vector
wmin = .350
wmax = .750
wVec = np.linspace(wmin,wmax,nW)

# Height vector
hmin = .150
hmax = .750
hVec = np.linspace(hmin,hmax,nH)

# Wavelength vector
lambdaMin = 1.4
lambdaMax = 1.7
lambdaVec = np.linspace(lambdaMin,lambdaMax,nLambda)

# Length vector
lMax = 150
nL   = 2
lVec = np.linspace(lMin,lMax,nL)

# Loss
alpha = 2 # dB/cm
alpha_ac = alpha / 4.34
alpha_ac_microns = alpha_ac * 1e-4

sc_y = 2  # supercell width (um)
sc_z = 2  # supercell height (um)
geometry_lattice = mp.Lattice(size=mp.Vector3(0, sc_y, sc_z))
resolution = 32  # pixels/um

ms = mpb.ModeSolver(
    geometry_lattice=geometry_lattice,
    resolution=resolution,
)

# Preallocate
neff = np.zeros([nW,nH,nLambda], dtype='float64')
S12  = np.zeros([nL,nW,nH,nLambda], dtype='complex128')

# ------------------------------------------------------------------------------ #
# Loop through experiment
# ------------------------------------------------------------------------------ #

for iW in range(nW):
    for iH in range(nH):
        for iLambda in range(nLambda):

            # get current frequency and wavelength
            currentLambda = lambdaVec[iLambda]
            omegaIn       = 1 / currentLambda

            # Update material parameters
            Si = mp.Medium(epsilon=epsSi(currentLambda))

            # Update geometry
            wCurrent = wVec[iW]
            hCurrent = hVec[iH]
            ms.geometry = [mp.Block(size=mp.Vector3(mp.inf, wCurrent, hCurrent), material=Si)]

            # Find k vector
            k = np.array(
                ms.find_k(
                    p              = mp.NO_PARITY,       # Polarization
                    omega          = omegaIn,            # Omega to find corresponding k
                    band_min       = 1,                  # Minimum band index
                    band_max       = 1,                  # Max band index
                    korig_and_kdir = mp.Vector3(1),      # K vector orientation
                    tol            = 1e-3,               # solver tolerance
                    kmag_guess     = omegaIn * 3.45,     # initial guess
                    kmag_min       = omegaIn * 0.1,      # Minimum
                    kmag_max       = omegaIn * 4))       # Maximum

            neff[iW,iH,iLambda]  = k / omegaIn
            S12[:,iW,iH,iLambda] = np.exp(-alpha_ac_microns * lVec) * np.exp(1j * 2 * np.pi * k * lVec)

hf = h5py.File('waveguide.h5', 'w')
hf.create_dataset('neff', data=neff, dtype='float64')
hf.create_dataset('S12', data=S12, dtype='complex128')
hf.close()
