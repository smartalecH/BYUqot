# ---------------------------------------------------------------------------- #
# sparam.py
# ---------------------------------------------------------------------------- #
# A simple script that illustrates how to extract s parameters from a simple
# waveguide using meep.
#
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

import meep as mp
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from matplotlib.widgets import Slider, Button, RadioButtons

# ---------------------------------------------------------------------------- #
# Important constants and materials
# ---------------------------------------------------------------------------- #

'''
eps0 = 1.3
susceptibilities = [mp.DrudeSusceptibility(frequency=5.28366673e-01, gamma=1.62316031, sigma=1.17034609e-03),
                    mp.LorentzianSusceptibility(frequency=4.19962269e+02, gamma=1.84154007e+03, sigma=7.44917949e+02)]
FR4 = mp.Medium(epsilon=eps0, E_susceptibilities=susceptibilities)
'''


Metal     = mp.Medium(epsilon=10)
Substrate = mp.Medium(epsilon=2)

resolution = 10                    # computational grid resolution

# ---------------------------------------------------------------------------- #
# Create Geometry
# ---------------------------------------------------------------------------- #

# Fundamental unit = 1 mm

# ------------------- Define substrate -------------------- #

substrateWidth     = 2.5
substrateHeight    = 10
substrateThickness = 0.25

SubstrateBlock = mp.Block(
    material = Substrate,
    center=mp.Vector3(0,0,0),
    size=mp.Vector3(substrateWidth,substrateHeight,substrateThickness)
)

# ------------------- Define square structures -------------------- #

metalThickness = 1/resolution;
boxZ = substrateThickness/2

linewidth = 0.2
outerBoxOuterWidth = 2.2
outerBoxInnerWidth = outerBoxOuterWidth - 2*linewidth

outerBoxMetal = mp.Block(
    material = Metal,
    center=mp.Vector3(0,0,boxZ),
    size=mp.Vector3(outerBoxOuterWidth,outerBoxOuterWidth,metalThickness)
)
outerBoxSub = mp.Block(
    material = Substrate,
    center=mp.Vector3(0,0,boxZ),
    size=mp.Vector3(outerBoxInnerWidth,outerBoxInnerWidth,metalThickness)
)

ringSeparation = 0.15
innerBoxOuterWith = outerBoxInnerWidth - 2*ringSeparation
innerBoxInnerWith = innerBoxOuterWith - 2*linewidth

innerBoxMetal = mp.Block(
    material = Metal,
    center=mp.Vector3(0,0,boxZ),
    size=mp.Vector3(innerBoxOuterWith,innerBoxOuterWith,metalThickness)
)
innerBoxSub = mp.Block(
    material = Substrate,
    center=mp.Vector3(0,0,boxZ),
    size=mp.Vector3(innerBoxInnerWith,innerBoxInnerWith,metalThickness)
)

ringGap = 0.3
ringTopY = -(outerBoxOuterWidth/2 - ringGap/2)
gapTop = mp.Block(
    material = Substrate,
    center=mp.Vector3(0,ringTopY,boxZ),
    size=mp.Vector3(ringGap,2*linewidth,metalThickness)
)

ringGap = 0.3
ringBottomY = (innerBoxOuterWith/2 - ringGap/2)
gapBottom=  mp.Block(
    material = Substrate,
    center=mp.Vector3(0,ringBottomY,boxZ),
    size=mp.Vector3(ringGap,2*linewidth,metalThickness)
)

# ------------------- Define wires -------------------- #

wireThickness = 0.14
wireLength = 10
wireZ = -substrateThickness/2
centerLine = mp.Block(
    material = Metal,
    center=mp.Vector3(0,0,wireZ),
    size=mp.Vector3(wireThickness,wireLength,metalThickness)
)

# ------------------- Consolidate geometries -------------------- #

geometry = [SubstrateBlock,outerBoxMetal,outerBoxSub,innerBoxMetal,innerBoxSub,gapTop,
            gapBottom,centerLine]

# ---------------------------------------------------------------------------- #
# Setup Simulation
# ---------------------------------------------------------------------------- #

# ------------------- Computational cell size -------------------- #
sx = 10   # computational grid width
sy = 2.5  # computaitonal grid height
sz = 2.5  # computaitonal grid thickness
cell_size = mp.Vector3(sx, sy, sz)  # computational grid

# ------------------- PML boundary layers -------------------- #
dpml = 1.0  # size of boundary layer
boundary_layers = [
    mp.PML(dpml,direction=mp.Y),
    mp.PML(dpml,direction=mp.Z) ,
    ]

# ------------------- Source -------------------- #
fLower = 1/75                   # lower frequency
fUpper = 1/15                   # upper frequency
fcen = np.mean([fLower,fUpper]) # center frequency
fwidth = fUpper - fLower        # frequency width

sources = [mp.Source(mp.GaussianSource(frequency=fcen,width=fwidth),
                     component=mp.Ez,
                     center=mp.Vector3(0,0))]

# ------------------- Simulation Block -------------------- #
sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=boundary_layers,
                    geometry=geometry,
                    k_point = mp.Vector3(1,0,0),
                    sources=sources)

# ---------------------------------------------------------------------------- #
# Extract S Parameters
# ---------------------------------------------------------------------------- #

nfreq = 4 # Number of frequency bins to compute

# ------------------- Port 1 Monitor -------------------- #

xm1 = -0.5*sx  # x-coordinate of monitor
mflux1 = sim.add_eigenmode(fcen, fwidth, nfreq,
                          mp.FluxRegion(center=mp.Vector3(xm1,0,0),
                                        size=mp.Vector3(0,sy-2*dpml,sz-2*dpml),
                                        direction=mp.X))

# ------------------- Port 2 Monitor -------------------- #

xm2 = 0.5*sx  # x-coordinate of monitor
mflux2 = sim.add_eigenmode(fcen, fwidth, nfreq,
                          mp.FluxRegion(center=mp.Vector3(xm2,0,0),
                                        size=mp.Vector3(0,sy-2*dpml,sz-2*dpml),
                                        direction=mp.X))

# ------------------- Run Simulation -------------------- #

sim.run(until=1)
#sim.run(until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(xm,0), 1e-10))

# ------------------- Extract S Params -------------------- #

bands = [1]  # indices of modes for which to compute expansion coefficients
alpha = sim.get_eigenmode_coefficients(mflux1, bands)
print(alpha)
'''
alpha0Plus  = alpha[2*0 + 0]  # coefficient of forward-traveling fundamental mode
alpha0Minus = alpha[2*0 + 1]  # coefficient of backward-traveling fundamental mode

print("refl:, {}, {:.8f}".format(Lt, abs(alpha0Minus)**2))
'''
# ---------------------------------------------------------------------------- #
# Plot and Compare
# ---------------------------------------------------------------------------- #

#eps_data = sim.get_array(center=mp.Vector3(), size=cell_size, component=mp.Dielectric)


'''
sio.savemat('eps_data.mat', mdict={'eps_data': eps_data})

celldims = eps_data.shape
numX = celldims[0]
numY = celldims[1]
numZ = celldims[2]

currentLayer = eps_data[:,:,0]

# generate figure
fig = plt.figure()
ax = plt.subplot(111)
fig.subplots_adjust(left=0.25, bottom=0.25)

# display image
l = ax.imshow(currentLayer.transpose(), interpolation='spline36', cmap='binary')

axcolor = 'lightgoldenrodyellow'
axfreq = fig.add_axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

sfreq = Slider(axfreq, 'Freq', 0, numZ-1, valinit=0, valstep=1)
print(numZ)

def update(val):
    currentLayer = eps_data[:,:,int(sfreq.val)]
    l.set_data(currentLayer.transpose())
    print(np.max(currentLayer))
    #draw()
    fig.canvas.draw()

sfreq.on_changed(update)

plt.figure(dpi=100)
plt.subplot(121)

maxEps = eps_data.sum(axis=0).sum(axis=1)
print(np.max(eps_data))
maxEps = maxEps[::-1].argsort()

print(maxEps)
frontFace = eps_data[:,:,69]
plt.imshow(frontFace.transpose(), interpolation='spline36', cmap='binary')
#plt.axis('off')

plt.subplot(122)
backFace = eps_data[:,:,56]
plt.imshow(backFace.transpose(), interpolation='spline36', cmap='binary')
#plt.axis('off')
plt.show()
'''
