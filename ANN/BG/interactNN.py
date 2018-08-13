from numpy import pi, sin
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from keras import backend as K
from keras.models import load_model
from numpy import matlib as npm
from numpy import linalg as LA
import pickle

# ---------------------------------------------------------------------------- #
# NN Function
# ---------------------------------------------------------------------------- #

lambda0   = 1.45         # Initial wavelength (microns)
lambdaEnd = 1.65         # Final wavelength (microns)
nLambda   = 1000         # Number of wavelength points
wavelengths = np.linspace(lambda0,lambdaEnd,nLambda)

with open('NN_bg_norm_1000.pkl', 'rb') as f:
    data = pickle.load(f)

datax = data[0]
magY = data[1]
aglY = data[2]


def r2Metric(y_true, y_pred):
    SS_res =  K.sum(K.square( y_true-y_pred ))
    SS_tot = K.sum(K.square( y_true - K.mean(y_true) ) )
    return ( 1 - SS_res/(SS_tot + K.epsilon()) )

# Load in the model
modelMag = load_model('NN_bg_mag_1000.h5',custom_objects={'r2Metric': r2Metric})
modelAgl = load_model('NN_bg_agl_1000.h5',custom_objects={'r2Metric': r2Metric})

def evalNN(a1,a2,NG,dw):
    a1Vec = np.squeeze(npm.repmat(a1,nLambda,1))
    a2Vec = np.squeeze(npm.repmat(a2,nLambda,1))
    NGVec = np.squeeze(npm.repmat(NG,nLambda,1))
    dwVec = np.squeeze(npm.repmat(dw,nLambda,1))

    INPUT = np.vstack((wavelengths,a1Vec,a2Vec,NGVec,dwVec)).T

    Mag   = magY.inverse_transform(np.reshape(modelMag.predict(datax.transform(INPUT)),(-1,1)))
    Agl   = aglY.inverse_transform(np.reshape(modelAgl.predict(datax.transform(INPUT)),(-1,1)))

    return Mag, Agl

# ---------------------------------------------------------------------------- #
# Plots
# ---------------------------------------------------------------------------- #

axis_color = 'lightgoldenrodyellow'

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

# Adjust the subplots region to leave some space for the sliders and buttons
fig.subplots_adjust(bottom=0.35)

a1_0 = 0.16
a2_0 = 0.16
NG_0 = 250
dw_0 = 0.05

mag0, vgrp0 = evalNN(a1_0, a2_0, NG_0, dw_0)

# Draw the initial plot
# The 'line' variable is used for modifying the line later
[line1] = ax1.plot(wavelengths, mag0, linewidth=2, color='red')
ax1.set_ylim([0, 1])

[line2] = ax2.plot(wavelengths, vgrp0, linewidth=2, color='red')
ax2.set_ylim([-1, 6])


# ---------------------------------------------------------------------------- #
# Slider bars
# ---------------------------------------------------------------------------- #

# Define an axes area and draw a slider in it
a1_slider_ax  = fig.add_axes([0.15, 0.25, 0.65, 0.03])
a1_slider = Slider(a1_slider_ax, 'a1', 0.15, 0.18, valinit=a1_0)

# Draw another slider
a2_slider_ax = fig.add_axes([0.15, 0.2, 0.65, 0.03])
a2_slider = Slider(a2_slider_ax, 'a2', 0.15, 0.18, valinit=a2_0)

# Draw another slider
NG_slider_ax = fig.add_axes([0.15, 0.15, 0.65, 0.03])
NG_slider = Slider(NG_slider_ax, 'NG', 100, 2000, valinit=NG_0)

# Draw another slider
dw_slider_ax = fig.add_axes([0.15, 0.1, 0.65, 0.03])
dw_slider = Slider(dw_slider_ax, 'dw', 0.01, 0.1, valinit=dw_0)

# ---------------------------------------------------------------------------- #
# Action calls
# ---------------------------------------------------------------------------- #

# Define an action for modifying the line when any slider's value changes
def sliders_on_changed(val):
    mag0, vgrp0 = evalNN(a1_slider.val, a2_slider.val, NG_slider.val, dw_slider.val)
    line1.set_ydata(mag0)
    line2.set_ydata(vgrp0)
    fig.canvas.draw_idle()

a1_slider.on_changed(sliders_on_changed)
a2_slider.on_changed(sliders_on_changed)
NG_slider.on_changed(sliders_on_changed)
dw_slider.on_changed(sliders_on_changed)


plt.show()
