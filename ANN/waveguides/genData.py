
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
# Define size of simulation array
# ------------------------------------------------------------------------------ #

# Extract the arguments passed
inputVec  = np.array(sys.argv[1:], dtype='float64')
wIter     = int(inputVec[0])
hIter     = int(inputVec[1])
numLambda = int(inputVec[2])

nW      = int(200)
nH      = int(50)
nLambda = numLambda

# Width vector
wmin = .350
wmax = 1.5
wVec = np.linspace(wmin,wmax,nW)

# Height vector
hmin = .150
hmax = .400
hVec = np.linspace(hmin,hmax,nH)

# Wavelength vector
lambdaMin = 1.4
lambdaMax = 1.7
lambdaVec = np.linspace(lambdaMin,lambdaMax,nLambda)

wCurrent = wVec[wIter]
hCurrent = hVec[hIter]

# ------------------------------------------------------------------------------ #
# Setup simulation domain
# ------------------------------------------------------------------------------ #

c = 299792458    # Speed of light in vacuum, m/s

# Silicon dispersion relation
eps0Si       = 7.98737492
epsLorentzSi = 3.68799143
omega0Si     = 3.93282466e15
epsSi        = lambda lam: eps0Si + epsLorentzSi*omega0Si ** 2 / (omega0Si ** 2 - (2*np.pi*c*1e6/lam)**2)

# Silicon Dioxide dispersion relation
eps0SiO2       = 2.119881
epsLorentzSiO2 = 49.43721
omega0SiO2     = 3.309238e13
epsSiO2        = lambda lam: eps0SiO2 + epsLorentzSiO2*omega0SiO2 ** 2 / (omega0SiO2 ** 2 - (2*np.pi*c*1e6/lam)**2)

sc_y = 2  # supercell width (um)
sc_z = 2  # supercell height (um)
geometry_lattice = mp.Lattice(size=mp.Vector3(0, sc_y, sc_z))
resolution = 32  # pixels/um

numModes = 2

ms = mpb.ModeSolver(
    geometry_lattice=geometry_lattice,
    resolution=resolution,
    ensure_periodicity=False
)

# Preallocate
neff_TE = np.zeros([nLambda,numModes], dtype='float64')
neff_TM = np.zeros([nLambda,numModes], dtype='float64')

# ------------------------------------------------------------------------------ #
# Loop through experiment
# ------------------------------------------------------------------------------ #

for iLambda in range(nLambda):

    # get current frequency and wavelength
    currentLambda = lambdaVec[iLambda]
    omegaIn       = 1 / currentLambda

    print('################################################')
    print(currentLambda)
    print('################################################')

    # Update material parameters
    Si   = mp.Medium(epsilon=epsSi(currentLambda))
    SiO2 = mp.Medium(epsilon=epsSiO2(currentLambda))

    # Update geometry
    ms.geometry = [mp.Block(size=mp.Vector3(mp.inf, wCurrent, hCurrent), material=Si)]
    ms.default_material = SiO2

    # Find TE modes vector
    kTE = np.array(
        ms.find_k(
            p              = mp.EVEN_Z,          # Polarization
            omega          = omegaIn,            # Omega to find corresponding k
            band_min       = 1,                  # Minimum band index
            band_max       = numModes,           # Max band index
            korig_and_kdir = mp.Vector3(1),      # K vector orientation
            tol            = 1e-3,               # solver tolerance
            kmag_guess     = omegaIn * 3.45,     # initial guess
            kmag_min       = omegaIn * 0.1,      # Minimum
            kmag_max       = omegaIn * 4))       # Maximum

    # Find TM modes
    kTM = np.array(
        ms.find_k(
            p              = mp.ODD_Z,           # Polarization
            omega          = omegaIn,            # Omega to find corresponding k
            band_min       = 1,                  # Minimum band index
            band_max       = numModes,           # Max band index
            korig_and_kdir = mp.Vector3(1),      # K vector orientation
            tol            = 1e-3,               # solver tolerance
            kmag_guess     = omegaIn * 3.45,     # initial guess
            kmag_min       = omegaIn * 0.1,      # Minimum
            kmag_max       = omegaIn * 4))       # Maximum
    # Store the effective indices for each mode
    neff_TE[iLambda,:]  = kTE / omegaIn
    neff_TM[iLambda,:]  = kTM / omegaIn

# ------------------------------------------------------------------------------ #
# Store results
# ------------------------------------------------------------------------------ #

filename = 'data_neff/wg_w-' + str(wIter) + '_h-' + str(hIter) + '.h5'

hf = h5py.File(filename, 'w')
hf.create_dataset('neff_TE', data=neff_TE, dtype='float64')
hf.create_dataset('neff_TM', data=neff_TM, dtype='float64')
hf.create_dataset('lambdaVec', data=lambdaVec, dtype='float64')
hf.create_dataset('wCurrent', data=wCurrent, dtype='float64')
hf.create_dataset('hCurrent', data=hCurrent, dtype='float64')
hf.close()

print('################################################')
print('Job finished')
print('################################################')
