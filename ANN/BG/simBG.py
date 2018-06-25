# ---------------------------------------------------------------------------- #
# simBG.py
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# Import libraries
# ---------------------------------------------------------------------------- #

import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize as op
from numpy import linalg as LA
import time

# ---------------------------------------------------------------------------- #
# Grating parameters
# ---------------------------------------------------------------------------- #

# Grating Parameters
Period    = 310e-9;       # Bragg period
NG        = 1000;          # Number of grating periods
L         = NG*Period;    # Grating length
width0    = 0.5;          # mean waveguide width
dwidth    = 0.005;        # +/-  waveguide width
width1    = width0 - dwidth;
width2    = width0 + dwidth;
loss_dBcm = 3;            # waveguide loss, dB/cm
loss      = np.log(10)*loss_dBcm/10*100;

# Simulation Parameters:
span      = 30e-9;        # Set the wavelength span for the simulation
Npoints   = 30000;

#500x220 oxide strip waveguide
neff_wavelength = lambda w: 2.4379 - 1.1193 * (w*1e6-1.554) - 0.0350 * (w*1e6-1.554) ** 2;
dneff_width     = lambda w: 10.4285*(w-0.5) ** 3 - 5.2487*(w-0.5) ** 2 + 1.6142*(w-0.5);

# Find Bragg wavelength using lambda_Bragg = Period * 2neff(lambda_bragg);
# Assume neff is for the average waveguide width.
f           = lambda lambdaIn: lambdaIn - Period*2* (neff_wavelength(lambdaIn)+(dneff_width(width2)+dneff_width(width1))/2);
wavelength0 = op.fsolve(f,1550e-9);
wavelengths = wavelength0 + np.linspace(-span/2, span/2, Npoints);
n1          = neff_wavelength(wavelengths)+dneff_width(width1);  # low index
n2          = neff_wavelength(wavelengths)+dneff_width(width2);  # high index

# ---------------------------------------------------------------------------- #
# Grating functions
# ---------------------------------------------------------------------------- #

def TMM_Grating_RT(wavelength, Period, NG, n1, n2, loss):
    # Calculate the R and T versus wavelength
    M = TMM_Grating_Matrix(wavelength, Period, NG, n1, n2, loss);
    S = np.zeros(M.shape,dtype='Complex128')
    S[0,0,:] = - M[1,0,:] / M[0,0,:]
    S[0,1,:] = 1 / M[0,0,:]
    S[0,0,:] = M[0,0,:] - M[0,1,:] * M[1,0,:] / M[0,0,:]
    S[0,0,:] = M[0,1,:] / M[0,0,:]

    q = wavelength.shape[0];
    T = np.abs(np.ones((q)) / np.squeeze(M[0,0,:])) ** 2;
    R = np.abs(np.squeeze(M[1,0,:]) / np.squeeze(M[0,0,:])) ** 2;
    return R, T, S

def TMM_Grating_Matrix(wavelength, Period, NG, n1, n2, loss):
    # Calculate the total transfer matrix of the gratings
    l      = Period/2;
    T_hw1  = TMM_HomoWG_Matrix(wavelength,l,n1,loss);
    T_is12 = TMM_IndexStep_Matrix(n1,n2);
    T_hw2  = TMM_HomoWG_Matrix(wavelength,l,n2,loss);
    T_is21 = TMM_IndexStep_Matrix(n2,n1);
    q      = wavelength.shape[0];
    Tp     = np.zeros((2,2,q),dtype='Complex128');
    T      = Tp;
    for i in range(wavelength.shape[0]):
    	Tp[:,:,i] = np.mat(T_hw2[:,:,i]) * np.mat(T_is21[:,:,i]) * np.mat(T_hw1[:,:,i]) * np.mat(T_is12[:,:,i])
    	T[:,:,i] = LA.matrix_power(np.mat(Tp[:,:,i]), NG); # 1st order uniform Bragg grating
    return T

def TMM_HomoWG_Matrix(wavelength,l,neff,loss):
    # Calculate the transfer matrix of a homogeneous waveguide.
    beta        = 2*np.pi*neff / wavelength-1j*loss/2; # Complex propagation constant
    T_hw        = np.zeros((2,2,neff.shape[0]),dtype='Complex128');
    T_hw[0,0,:] = np.exp(1j*beta*l);
    T_hw[1,1,:] = np.exp(-1j*beta*l);
    return T_hw

def TMM_IndexStep_Matrix(n1,n2):
    # Calculate the transfer matrix for a index step from n1 to n2.
    T_is        = np.zeros((2,2,n1.shape[0]),dtype='Complex128');
    a           = (n1+n2) / (2*np.sqrt(n1*n2));
    b           = (n1-n2) / (2*np.sqrt(n1*n2));
    T_is[0,0,:] = a;
    T_is[0,1,:] = b;
    T_is[1,0,:] = b;
    T_is[1,1,:] = a;
    return T_is

# ---------------------------------------------------------------------------- #
# Grating functions
# ---------------------------------------------------------------------------- #
T_1 = time.clock()
R,T,S = TMM_Grating_RT(wavelengths, Period, NG, n1, n2, loss);
R1 = abs(S[1,1,:]) ** 2
T_2 = time.clock()

print(T_2 - T_1)

plt.figure();
plt.plot(wavelengths*1e6,R, linewidth=2.0);
plt.plot(wavelengths*1e6,R1, '--',linewidth=2.0);
plt.xlabel('Wavelength $\mu m$')
plt.ylabel('Response');

plt.show()
