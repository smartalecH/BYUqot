# ------------------------------------------------------------------ #
#       SiP_Materials.py
# ------------------------------------------------------------------ #
#
# Various materials used in silicon photonics
#
# ------------------------------------------------------------------ #
#       VERSION HISTORY
# ------------------------------------------------------------------ #
# 18 May 2018 - AMH - Initialization
#
# ------------------------------------------------------------------ #
#       To - DO
# ------------------------------------------------------------------ #
#
# ------------------------------------------------------------------ #
# Import packages
# ------------------------------------------------------------------ #

import meep as mp
import numpy as np

# ------------------------------------------------------------------ #
# Relevant Parameters
# ------------------------------------------------------------------ #

# default unit length is 1 um
um_scale = 1.0

# conversion factor for eV to 1/um [=1/hc]
eV_um_scale = um_scale/1.23984193

# Speed of light in vacuum, m/s
c = 299792458

# ------------------------------------------------------------------ #
# Silicon (Si)
# ------------------------------------------------------------------ #
# silicon (Si) from Palik's Handbook & Lukas's book
# wavelength range: 1.15 - 1.8 um

# Reference:
# H. H. Li. Refractive index of silicon and germanium and its wavelength and
# temperature derivatives, J. Phys. Chem. Ref. Data 9, 561-658 (1993)

Si_range = mp.FreqRange(min=1/1.8, max=1/1.15)

eps = 7.9874
eps_lorentz = 3.6880
omega0 = 3.9328e15
delta0 = 0

Si_frq1 = omega0 / (2 * np.pi * c) * 1e-6
Si_gam1 = 0
Si_sig1 = eps_lorentz

Si_susc = [ mp.LorentzianSusceptibility(frequency=Si_frq1, gamma=Si_gam1, sigma=Si_sig1) ]

Si = mp.Medium(epsilon=eps, E_susceptibilities=Si_susc, valid_freq_range=Si_range)

# ------------------------------------------------------------------ #
# Silicon Dioxide (SiO2)
# ------------------------------------------------------------------ #
# silicon dioxide (SiO2) from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.25 - 1.77 um

SiO2_range = mp.FreqRange(min=1/1.77, max=1/0.25)

SiO2_frq1 = 1/(0.103320160833333*um_scale)
SiO2_gam1 = 1/(12.3984193000000*um_scale)
SiO2_sig1 = 1.12

SiO2_susc = [ mp.LorentzianSusceptibility(frequency=SiO2_frq1, gamma=SiO2_gam1, sigma=SiO2_sig1) ]

SiO2 = mp.Medium(epsilon=1.0, E_susceptibilities=SiO2_susc, valid_freq_range=SiO2_range)
# ------------------------------------------------------------------ #
# Silicon Nitride non-stoichiometric
# ------------------------------------------------------------------ #

# silicon nitride (SiN), non-stoichiometric, from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.21 - 2.07 um

SiN_range = mp.FreqRange(min=1/2.07, max=1/0.21)

SiN_frq1 = 1/(0.190891752117013*um_scale)
SiN_gam1 = 1/(3.11518072864322*um_scale)
SiN_sig1 = 1.2650

SiN_susc = [ mp.LorentzianSusceptibility(frequency=SiN_frq1, gamma=SiN_gam1, sigma=SiN_sig1) ]

SiN = mp.Medium(epsilon=2.320, E_susceptibilities=SiN_susc, valid_freq_range=SiN_range)

# ------------------------------------------------------------------ #
# Silicon Nitride stoichiometric
# ------------------------------------------------------------------ #

# silicon nitride (Si3N4), stoichiometric, from Horiba Technical Note 08: Lorentz Dispersion Model
# ref: http://www.horiba.com/fileadmin/uploads/Scientific/Downloads/OpticalSchool_CN/TN/ellipsometer/Lorentz_Dispersion_Model.pdf
# wavelength range: 0.23 - 0.83 um

Si3N4_range = mp.FreqRange(min=1/0.83, max=1/0.23)

Si3N4_frq1 = 1/(0.389153148148148*um_scale)
Si3N4_gam1 = 1/(0.693811936205932*um_scale)
Si3N4_sig1 = 4.377

Si3N4_susc = [ mp.LorentzianSusceptibility(frequency=Si3N4_frq1, gamma=Si3N4_gam1, sigma=Si3N4_sig1) ]

Si3N4 = mp.Medium(epsilon=1.0, E_susceptibilities=Si3N4_susc, valid_freq_range=Si3N4_range)
