import meep as mp
from meep import mpb
import numpy as np

import matplotlib.pyplot as plt

# Compute modes of a rectangular Si strip waveguide on top of oxide.
# Note that you should only pay attention, here, to the guided modes,
# which are the modes whose frequency falls under the light line --
# that is, frequency < beta / 1.45, where 1.45 is the SiO2 index.

# Since there's no special lengthscale here, I'll just
# use microns.  In general, if you use units of x, the frequencies
# output are equivalent to x/lambda# so here, the freqeuncies will be
# output as um/lambda, e.g. 1.5um would correspond to the frequency
# 1/1.5 = 0.6667.

w = 0.5  # Si width (um)
h = 0.22  # Si height (um)

Si = mp.Medium(index=3.4757)
SiO2 = mp.Medium(index=1.444)

# Define the computational cell.  We'll make x the propagation direction.
# the other cell sizes should be big enough so that the boundaries are
# far away from the mode field.
sc_y = 2  # supercell width (um)
sc_z = 2  # supercell height (um)

geometry_lattice = mp.Lattice(size=mp.Vector3(0, sc_y, sc_z))
# define the 2d blocks for the strip and substrate
width     = 0.5
height    = 0.22
thickness = 1e20
waveguide = mp.Block(
    size=mp.Vector3(thickness,width,height),
    center=mp.Vector3(0,0,0),
    material=Si)
geometry = [waveguide]

# The k (i.e. beta, i.e. propagation constant) points to look at, in
# units of 2*pi/um.  We'll look at num_k points from k_min to k_max.
num_k = 9
k_min = 0.1
k_max = 3.0
k_points = mp.interpolate(num_k, [mp.Vector3(k_min), mp.Vector3(k_max)])

k_points = [mp.Vector3(1/1.55)]
resolution = 64  # pixels/um

# Increase this to see more modes.  (The guided ones are the ones below the
# light line, i.e. those with frequencies < kmag / 1.45, where kmag
# is the corresponding column in the output if you grep for "freqs:".)


filename_prefix = 'strip-'  # use this prefix for output files

ms = mpb.ModeSolver(
    geometry_lattice=geometry_lattice,
    geometry=geometry,
    resolution=resolution,
    filename_prefix=filename_prefix,
    default_material=SiO2,
)

# compute num_bands lowest frequencies as a function of k. Also display
# "parities", i.e. whether the mode is symmetric or anti_symmetric
# through the y=0 and z=0 planes.
#ms.run()

# Output the x component of the Poynting vector for num_bands bands at omega
num_bands = 1
efields = []
def addField(tr_ms, band):
    efields.append(tr_ms.get_efield(band))

lam = 1.55
omega = 1 / 1.55  # frequency corresponding to 1.55um

foundk = ms.find_k(mp.EVEN_Z,
                    omega,
                    1,
                    num_bands,
                    mp.Vector3(1),
                    1e-3,
                    omega * 3.45,
                    omega * 0.1,
                    omega * 4,
                    addField)
foundk = np.array(foundk)
n = foundk * lam
print(n)
first_band = 1
#temp = ms.get_poynting(first_band,num_bands)
efield = ms.get_efield(1)
#field2 = np.squeeze(efields[2])
#print(field1.shape)


temp = efield[:,:,0,:]
#temp = field2
Ex = temp[:,:,0]
Ey = temp[:,:,1]
Ez = temp[:,:,2]

I = (abs(Ex)**2 + abs(Ey)**2 + abs(Ez)**2)
Ix = abs(Ex)**2
Iy = abs(Ey)**2
Iz = abs(Ez)**2


eps_data = ms.get_epsilon()

plt.figure(dpi=100)

plt.subplot(221)
plt.imshow(np.abs(eps_data.transpose()),interpolation='spline36',cmap='binary')
plt.imshow(np.abs(Ix.transpose()),interpolation='spline36',alpha=0.7)
plt.title('$I_x$')
plt.colorbar()

plt.subplot(222)
plt.imshow(np.abs(eps_data.transpose()),interpolation='spline36',cmap='binary')
plt.imshow(np.abs(Iy.transpose()),interpolation='spline36',alpha=0.7)
plt.title('$I_y$')
plt.colorbar()

plt.subplot(223)
plt.imshow(np.abs(eps_data.transpose()),interpolation='spline36',cmap='binary')
plt.imshow(np.abs(Iz.transpose()),interpolation='spline36',alpha=0.7)
plt.title('$I_z$')
plt.colorbar()

plt.subplot(224)
plt.imshow(np.abs(eps_data.transpose()),interpolation='spline36',cmap='binary')
plt.imshow(np.abs(I.transpose()),interpolation='spline36',alpha=0.7)
plt.title('$I$')
plt.colorbar()
plt.show()


###########################################################################

# Above, we outputted the dispersion relation: frequency (omega) as a
# function of wavevector kx (beta).  Alternatively, you can compute
# beta for a given omega -- for example, you might want to find the
# modes and wavevectors at a fixed wavelength of 1.55 microns.  You
# can do that using the find_k function:



'''

'''
