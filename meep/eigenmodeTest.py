import meep as mp
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Example file illustrating an eigenmode source, generating a waveguide mode
# (requires recent MPB version to be installed before Meep is compiled)



# ---------------------------------------------------------------------------- #
# Setup Constants
# ---------------------------------------------------------------------------- #

# Speed of light in vacuum, m/s
c = 299792458

# default unit length is 1 um
um_scale = 1.0

# import silicon
Si   = mp.Medium(index=3.4757)#SiP_Materials.Si

# import SiO2
SiO2 = mp.Medium(index=1.4440)#SiP_Materials.SiO2

# ---------------------------------------------------------------------------- #
# Draw geometries
# ---------------------------------------------------------------------------- #

# Define slab waveguide dimensions (ignore 3rd dimension)
width     = 0.5
height    = 0.22
thickness = 1e20
waveguide = mp.Block(
    size=mp.Vector3(thickness,width,height),
    center=mp.Vector3(0,0,0),
    material=Si)
geometry = [waveguide]

# ---------------------------------------------------------------------------- #
# Setup simulation
# ---------------------------------------------------------------------------- #

# Specify size of bondary layers
boundaryCondWidth = 2.0
pml_layers = [mp.PML(boundaryCondWidth)]

# Specify cell sizes
cellWidth = 2*width + 2*boundaryCondWidth
cellHeight = 2*height + 2*boundaryCondWidth
cell = mp.Vector3(0,cellWidth,cellHeight)


# create a transparent source that excites a right-going waveguide mode
modeNumber = 2
sources = [
    mp.EigenModeSource(src=mp.ContinuousSource(1/1.55), direction=mp.X,
                       center=mp.Vector3(),eig_band=modeNumber)
]


force_complex_fields = True  # so we can get time-average flux

resolution = 20

sim = mp.Simulation(
    cell_size=cell,
    geometry=geometry,
    sources=sources,
    boundary_layers=pml_layers,
    force_complex_fields=force_complex_fields,
    resolution=resolution
)

sim.run(
    mp.at_beginning(mp.output_epsilon),
    mp.at_end(mp.output_png(mp.Ex, "-a yarg -A $EPS -S3 -Zc dkbluered", rm_h5=False)),
    until=200
)

flux1 = sim.flux_in_box(mp.X, mp.Volume(center=mp.Vector3(-2.0), size=mp.Vector3(1.8, 6)))
flux2 = sim.flux_in_box(mp.X, mp.Volume(center=mp.Vector3(2.0), size=mp.Vector3(1.8, 6)))

# averaged over y region of width 1.8
print("left-going flux = {}".format(flux1 / -1.8))

# averaged over y regison of width 1.8
print("right-going flux = {}".format(flux2 / 1.8))

eps_data = sim.get_array(center=mp.Vector3(),size=cell,component = mp.Dielectric)
ex_data = sim.get_array(center=mp.Vector3(),size=cell,component = mp.Ex)
ey_data = sim.get_array(center=mp.Vector3(),size=cell,component = mp.Ey)
ez_data = sim.get_array(center=mp.Vector3(),size=cell,component = mp.Ez)

I = (abs(ex_data)**2 + abs(ey_data)**2 + abs(ez_data)**2)
Ix = abs(ex_data)**2
Iy = abs(ey_data)**2
Iz = abs(ez_data)**2

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
