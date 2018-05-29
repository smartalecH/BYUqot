# ---------------------------------------------------------------------------- #
# fitMaterials.py
# ---------------------------------------------------------------------------- #
# TODO:
#
# ---------------------------------------------------------------------------- #
# VERSION
# ---------------------------------------------------------------------------- #
# 29 May 2018 - AMH - Initialization
#
# ---------------------------------------------------------------------------- #
# Import Libraries
# ---------------------------------------------------------------------------- #

import numpy as np
import numpy.matlib as npm
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import leastsq

# ---------------------------------------------------------------------------- #
# Generate a Lorentzian function
# ---------------------------------------------------------------------------- #

def genLorentzian(f,params):

    # Check if input parameters are good
    numParams = len(params)
    if numParams < 7:
        print(numParams)
        raise ValueError("numParams does not have enough parameters")
    elif (numParams - 4) % 3 != 0:
        raise ValueError("numParams either has too many or too little parameters")

    # Initialize lorentzian with epsilon and drude term
    eps0 = params[0]
    f0 = params[1]
    sigma0 = params[2]
    gamma0 = params[3]
    num = f0**2 * sigma0
    den = f**2 - 1j * f * gamma0 / (2 * np.pi)
    lorentzian = eps0 + num / den

    # Iterate through lorentzian susceptibilities
    numSuscept = int((numParams - 1) / 3)
    for paramIn in range(0,numSuscept):
        f0 = params[3*paramIn + 1]
        sigma = params[3*paramIn + 2]
        gamma = params[3*paramIn + 3]

        num = f0**2 * sigma
        den = f0**2 - f**2 - 1j * f * gamma / (2 * np.pi)
        lorentzian = lorentzian + num / den

    return lorentzian

# ---------------------------------------------------------------------------- #
# Generate a Lorentzian cost function
# ---------------------------------------------------------------------------- #

def lorentzianCost(f,realL,imagL,params):
    lorentzian = genLorentzian(f,params)
    realError = np.abs(np.real(lorentzian) - realL) ** 2
    imagError = np.abs(np.imag(lorentzian) - imagL) ** 2

    totalError = np.sum(realError + imagError)

    modError = lorentzian - (realL + 1j*imagL)
    tempError = modError.real **2 + modError.imag **2
    return (tempError)

# ---------------------------------------------------------------------------- #
# Fit data to a Lorentzian
# ---------------------------------------------------------------------------- #

def fitLorentzian(f,realL,imagL,params):

    cost = lambda x: lorentzianCost(f,realL,imagL,x)
    numSuscept = 1
    x0 = np.squeeze(np.ones(4 + 3*numSuscept))

    #res = leastsq(cost,x0)
    lb = np.squeeze(np.zeros(4 + 3*numSuscept))
    ub = np.squeeze(npm.repmat(1e6,1 + 3*numSuscept,1))
    bounds = np.transpose(np.array([lb,ub]))
    res = minimize(cost, x0, bounds = bounds, method='trust-krylov')
    return res.x


# ---------------------------------------------------------------------------- #
# Test
# ---------------------------------------------------------------------------- #

N = 100
flower = 1/75
fupper = 1/15
f = np.linspace(flower,fupper,N)
params = [1, 2,1.5,0]

realL = np.squeeze(npm.repmat(4.399999,N,1))
imagL = np.squeeze(npm.repmat(0.0879998,N,1))
res = fitLorentzian(f,realL,imagL,params)

plt.figure()
plt.subplot(211)
plt.plot(f,realL)
plt.plot(f,np.real(genLorentzian(f,res)))

plt.subplot(212)
plt.plot(f,imagL)
plt.plot(f,np.imag(genLorentzian(f,res)))

plt.show()
