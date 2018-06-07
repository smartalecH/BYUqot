# ---------------------------------------------------------------------------- #
# materialSweep.py
# ---------------------------------------------------------------------------- #
#
# Investigated Materials:
# * NaCl
# * Sucrose
# * Ethylene Glycol
# * Glycerol
# * DMSO
#
# TODO:
#
# ---------------------------------------------------------------------------- #
# VERSION HISTORY
# ---------------------------------------------------------------------------- #
# 20 MAY 2018 - AMH - Initialization
# 07 JUN 2018 - AMH - Added Broadband Sweeps
#
# ---------------------------------------------------------------------------- #
# Import Libraries
# ---------------------------------------------------------------------------- #

import numpy as np
import h5py
import matplotlib.pyplot as plt
import scipy.io as sio

# ---------------------------------------------------------------------------- #
# Load, Fit, and Plot Effective Index Simulation Data
# ---------------------------------------------------------------------------- #

mat_contents = sio.loadmat('neffExperiment_simple.mat')

neff_1         = mat_contents['neff_1']
neff_2         = mat_contents['neff_2']
waveguideBank = np.squeeze(np.array(mat_contents['waveguideBank']))
claddingBank   = mat_contents['claddingBank']

# Setup polynomial fitting
nPts = 100
nCladding = np.linspace(1,3,nPts)
polyOrder = 4

polycoeffs = np.zeros([int(len(waveguideBank)/2),(polyOrder)])

# Setup figure
plt.figure(1)
legendLabel = []

linecolors = ['b', 'g', 'r', 'c', 'm']

# Iterate through half of waveguide widths
for n in range(1, len(waveguideBank), 2):
    # Perform polynomial fitting
    z = np.polyfit(np.squeeze(claddingBank), np.squeeze(neff_1[n,:]), polyOrder)
    p = np.poly1d(z)
    nFit = p(nCladding)

    # Organize legend labels
    label = '%i nm' %  ((waveguideBank[n])*1.0e9)

    # Plot data
    plt.plot(claddingBank,neff_1[n,:],marker='o', color=linecolors[int(n/2)],linestyle='None',label=label)
    plt.plot(nCladding,nFit, color=linecolors[int(n/2)],label=None)

# Clean up plots
plt.legend()
plt.grid(True,color='0.85')
plt.xlabel('Index of Cladding')
plt.ylabel('Effective Index')
plt.savefig('effectiveIndex.png')

# ---------------------------------------------------------------------------- #
# Load Material Index Data
# ---------------------------------------------------------------------------- #
# Let's explore the empirically determined index evolution of certain
# solvents when dissolved in DI water.
#
# All data is pulled from the following reference:
# J. E. Saunders, C. Sanders, H. Chen, and H.-P. Loock,
# “Refractive indices of common solvents and solutions at 1550 nm,”
# Applied Optics, vol. 55, no. 4, p. 947, Feb. 2016.

# Number of molarity/concentration points to plot
NW = 100

# Begin new figure
plt.figure(2)

# ------------------- NaCl -------------------------------

WmaxNaCl = 0.25      # Maximum solubility as a weight fraction
WNaCl = np.linspace(0,WmaxNaCl,NW)
# Polynomial fit coefficients
A = -0.0800
B = 0.0740
C = 0.1620
D = 1.3162

coef  = np.array([A, B, C, D])
pNaCl = np.poly1d(coef)
yNaCl = pNaCl(WNaCl)

plt.plot(WNaCl,yNaCl,label='NaCl')

# -------------------  Sucrose ---------------------------

WmaxSucrose = 0.80      # Maximum solubility as a weight fraction
WSucrose = np.linspace(0,WmaxSucrose,NW)
# Polynomial fit coefficients
A = 0.0330
B = 0.0456
C = 0.1416
D = 1.3166

coef  = np.array([A, B, C, D])
pSucrose = np.poly1d(coef)
ySucrose = pSucrose(WSucrose)

plt.plot(WSucrose,ySucrose,label='Sucrose')

# -------------------  Ethylene Glycol -------------------

WmaxEG= 1      # Maximum solubility as a weight fraction
WEG = np.linspace(0,WmaxEG,NW)
# Polynomial fit coefficients
A = -0.0223
B = 0.0321
C = 0.0910
D = 1.3166

coef  = np.array([A, B, C, D])
pEG = np.poly1d(coef)
yEG = pEG(WSucrose)

plt.plot(WEG,yEG,label='Ethylene Glycol')

# -------------------  Glycerol --------------------------

WmaxGlycerol= 1      # Maximum solubility as a weight fraction
WGlycerol = np.linspace(0,WmaxGlycerol,NW)
# Polynomial fit coefficients
A = -0.0216
B = 0.0512
C = 0.1110
D = 1.3165

coef  = np.array([A, B, C, D])
pGlycerol = np.poly1d(coef)
yGlycerol = pEG(WGlycerol)

plt.plot(WGlycerol,yGlycerol,label='Glycerol')

# ------------------- DMSO -------------------------------

WmaxDMSO= 1      # Maximum solubility as a weight fraction
WDMSO = np.linspace(0,WmaxDMSO,NW)
# Polynomial fit coefficients
A = -0.0730
B = 0.0998
C = 0.1189
D = 1.3172

coef  = np.array([A, B, C, D])
pDMSO = np.poly1d(coef)
yDMSO = pDMSO(WDMSO)

plt.plot(WDMSO,yDMSO,label='DMSO')


# ------------------- Clean up plots --------------------

plt.legend()
plt.grid(True,color='0.85')
plt.xlabel('Mass Fraction, W')
plt.ylabel('Refractive Index')
plt.savefig('materialIndex.png')

# ---------------------------------------------------------------------------- #
# Load Material Index Data
# ---------------------------------------------------------------------------- #
# Let's see how these different materials can modulate the effective index
# of the waveguide as a function of molar mass concentration. We'll stick with
# a 500 nm wide waveguide for simplicity

# Perform polynomial fitting
n = 4
z = np.polyfit(np.squeeze(claddingBank), np.squeeze(neff_1[n,:]), polyOrder)
pw500 = np.poly1d(z)

plt.figure(3)

# ------------------- NaCl -------------------------------
neffNaCl = pw500(yNaCl)
plt.plot(WNaCl,neffNaCl,label='NaCl')
# -------------------  Sucrose ---------------------------
neffSucrose = pw500(ySucrose)
plt.plot(WSucrose,neffSucrose,label='Sucrose')
# -------------------  Ethylene Glycol -------------------
neffEG = pw500(yEG)
plt.plot(WEG,neffEG,label='Ethylene Glycol')
# -------------------  Glycerol --------------------------
neffGlycerol = pw500(yGlycerol)
plt.plot(WGlycerol,neffGlycerol,label='Glycerol')
# ------------------- DMSO -------------------------------
neffDMSO = pw500(yDMSO)
plt.plot(WDMSO,neffDMSO,label='DMSO')

plt.legend()
plt.grid(True,color='0.85')
plt.xlabel('Mass Fraction, W')
plt.ylabel('Effective Index')

plt.savefig('materialEffectiveIndex.png')

# ---------------------------------------------------------------------------- #
# MZI Results
# ---------------------------------------------------------------------------- #
# We can now evaluate how a Mach-Zehnder Interfoeromter (MZI) will interact with
# the above relationships. We can assume that the modulating channels envelope
# the entire device, so that both arms have the same new effective index.

# Begin new figure
#fig, ax = plt.subplots(1, 1, sharex='col', sharey='row')
plt.figure()
lam = 1550e-9
L1 = 250e-6
deltaL = [100e-6]

for k in range(0,len(deltaL)):
    neff = 1
    L2 = L1 + deltaL[k]
    I = lambda Beta: 0.5 * (1+np.cos(Beta * L1 - Beta * L2))

    row = int(k/2)
    col = k % 2

    # ------------------- NaCl -------------------------------
    BetaNaCl = 2 * np.pi * neffNaCl/lam
    INaCl = I(BetaNaCl)
    plt.plot(WNaCl,INaCl,label='NaCl')
    # -------------------  Sucrose ---------------------------
    BetaSucrose = 2 * np.pi * neffSucrose/lam
    ISucrose = I(BetaNaCl)
    plt.plot(WSucrose,ISucrose,label='Sucrose')
    # -------------------  Ethylene Glycol -------------------
    BetaEG = 2 * np.pi * neffEG/lam
    IEG = I(BetaEG)
    plt.plot(WEG,IEG,label='Ethylene Glycol')
    # -------------------  Glycerol --------------------------
    BetaGlycerol = 2 * np.pi * neffGlycerol/lam
    IGlycerol = I(BetaGlycerol)
    plt.plot(WGlycerol,IGlycerol,label='Glycerol')
    # ------------------- DMSO -------------------------------
    BetaDMSO = 2 * np.pi * neffDMSO/lam
    IDMSO = I(BetaDMSO)
    plt.plot(WDMSO,IDMSO,label='DMSO')

    #plt.title('$\Delta L =$ %i $\mu m$' %  ((deltaL[k])*1.0e6))
    plt.legend( ('NaCl','Sucrose','Ethylene Glycol','Glycerol','DMSO'),bbox_to_anchor=(0., 1.02, 1., .102), loc='upper center',
               ncol=3, mode="expand", borderaxespad=0 )
    plt.xlabel('Mass Fraction, W')
    plt.ylabel('Power')
#fig.text(0.5, 0.04, 'Mass Fraction, W', ha='center')
#fig.text(0.04, 0.5, 'Power', va='center', rotation='vertical')
plt.grid(True)
plt.savefig('materialIndex.png')

# ---------------------------------------------------------------------------- #
# Broadband results for MZI
# ---------------------------------------------------------------------------- #
# Now let's look at a braodband spectrum for various MZI path length differences
# and various concentrations of each solute

neff = 1
lambda0 = 1.55e-6
deltaLambda = 15e-9
lambdaStart = lambda0 - deltaLambda/2
lambdaEnd   = lambda0 + deltaLambda/2
numLambda   = 250

lambdaVec = np.linspace(lambdaStart,lambdaEnd,numLambda)

deltaL = [100e-6]

# ------------------- NaCl -------------------------------
#fig, ax = plt.subplots(2, 2, sharex='col', sharey='row')
lgndEntries = []
plt.figure()
for k in range(0,len(deltaL)):
    for kc in range(0,NW,int(NW/4)):
        L2 = L1 + deltaL[k]
        I = lambda Beta: 0.5 * (1+np.cos(Beta * L1 - Beta * L2))

        row = int(k/2)
        col = k % 2
        BetaNaCl = 2 * np.pi * neffNaCl[kc]/lambdaVec
        INaCl = I(BetaNaCl)
        concentration = "W = %.2f" % WNaCl[kc]
        plt.title('NaCl $\Delta L =$ %i $\mu m$' %  ((deltaL[k])*1.0e6))
        plt.plot(lambdaVec*1e9,INaCl,label=concentration)
        plt.xlabel('Wavelength, $nm$')
        plt.ylabel('Power')
        #fig.text(0.5, 0.04, 'Wavelength, $nm$', ha='center')
        #fig.text(0.04, 0.5, 'Power', va='center', rotation='vertical')
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='upper center',
               ncol=3, mode="expand", borderaxespad=0 )
    plt.grid(True)
    plt.savefig('NaCl.png')
# -------------------  Sucrose ---------------------------
#fig, ax = plt.subplots(2, 2, sharex='col', sharey='row')
plt.figure()
for k in range(0,len(deltaL)):
    for kc in range(0,NW,int(NW/4)):
        L2 = L1 + deltaL[k]
        I = lambda Beta: 0.5 * (1+np.cos(Beta * L1 - Beta * L2))

        row = int(k/2)
        col = k % 2
        BetaSucrose = 2 * np.pi * neffSucrose[kc]/lambdaVec
        ISucrose = I(BetaSucrose)
        concentration = "W = %.2f" % WSucrose[kc]
        plt.plot(lambdaVec*1e9,ISucrose,label=concentration)
        plt.title('Sucrose $\Delta L =$ %i $\mu m$' %  ((deltaL[k])*1.0e6))
        plt.xlabel('Wavelength, $nm$')
        plt.ylabel('Power')
        #fig.text(0.5, 0.04, 'Wavelength, $nm$', ha='center')
        #fig.text(0.04, 0.5, 'Power', va='center', rotation='vertical')
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='upper center',
                   ncol=3, mode="expand", borderaxespad=0 )
        plt.grid(True)
        plt.savefig('Sucrose.png')

# -------------------  Ethylene Glycol -------------------
plt.figure()
for k in range(0,len(deltaL)):
    for kc in range(0,NW,int(NW/4)):
        L2 = L1 + deltaL[k]
        I = lambda Beta: 0.5 * (1+np.cos(Beta * L1 - Beta * L2))

        row = int(k/2)
        col = k % 2
        BetaEG = 2 * np.pi * neffEG[kc]/lambdaVec
        IEG = I(BetaEG)
        concentration = "W = %.2f" % WEG[kc]
        plt.plot(lambdaVec*1e9,IEG,label=concentration)
        plt.title('Ethylene Glycol $\Delta L =$ %i $\mu m$' %  ((deltaL[k])*1.0e6))
        plt.xlabel('Wavelength, $nm$')
        plt.ylabel('Power')
        #fig.text(0.5, 0.04, 'Wavelength, $nm$', ha='center')
        #fig.text(0.04, 0.5, 'Power', va='center', rotation='vertical')
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='upper center',
                   ncol=3, mode="expand", borderaxespad=0 )
        plt.grid(True)
        plt.savefig('EG.png')

# -------------------  Glycerol --------------------------
plt.figure()
for k in range(0,len(deltaL)):
    for kc in range(0,NW,int(NW/4)):
        L2 = L1 + deltaL[k]
        I = lambda Beta: 0.5 * (1+np.cos(Beta * L1 - Beta * L2))

        row = int(k/2)
        col = k % 2
        BetaGlycerol = 2 * np.pi * neffGlycerol[kc]/lambdaVec
        IGlycerol = I(BetaGlycerol)
        concentration = "W = %.2f" % WGlycerol[kc]
        plt.plot(lambdaVec*1e9,IGlycerol,label='Glycerol')
        plt.title('Glycerol $\Delta L =$ %i $\mu m$' %  ((deltaL[k])*1.0e6))
        plt.xlabel('Wavelength, $nm$')
        plt.ylabel('Power')
        #fig.text(0.5, 0.04, 'Wavelength, $nm$', ha='center')
        #fig.text(0.04, 0.5, 'Power', va='center', rotation='vertical')
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='upper center',
                   ncol=3, mode="expand", borderaxespad=0 )
        plt.grid(True)
        plt.savefig('glycerol.png')
# ------------------- DMSO -------------------------------
plt.figure()
for k in range(0,len(deltaL)):
    for kc in range(0,NW,int(NW/4)):
        L2 = L1 + deltaL[k]
        I = lambda Beta: 0.5 * (1+np.cos(Beta * L1 - Beta * L2))

        row = int(k/2)
        col = k % 2
        BetaDMSO = 2 * np.pi * neffDMSO[kc]/lambdaVec
        IDMSO = I(BetaDMSO)
        concentration = "W = %.2f" % WDMSO[kc]
        plt.plot(lambdaVec*1e9,IDMSO,label=concentration)
        plt.title('DMSO $\Delta L =$ %i $\mu m$' %  ((deltaL[k])*1.0e6))
        plt.xlabel('Wavelength, $nm$')
        plt.ylabel('Power')
        #fig.text(0.5, 0.04, 'Wavelength, $nm$', ha='center')
        #fig.text(0.04, 0.5, 'Power', va='center', rotation='vertical')
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='upper center',
                   ncol=3, mode="expand", borderaxespad=0 )
        plt.grid(True)
        plt.savefig('DMSO.png')

# ---------------------------------------------------------------------------- #
# Plotting Cleanup
# ---------------------------------------------------------------------------- #

plt.show()
