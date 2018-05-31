import meep as mp
import math

resolution = 20        # pixels/um

nglass = 1.5
glass = mp.Medium(index=nglass)

dsub = 3               # substrate thickness
dair = 3               # air padding between grating and pml
gp = 10.0              # grating period
gh = 0.5               # grating height
gdc = 0.5              # grating duty cycle

lambda_min = 0.4       # vacuum wavelength
lambda_max = 0.8       # vacuum wavelength
fmin = 1/lambda_max
fmax = 1/lambda_min
fcen = 0.5*(fmin+fmax)
df = fmax-fmin

dpml = lambda_max      # PML thickness
sx = dpml+dsub+gh+dair+dpml
sy = gp

cell_size = mp.Vector3(sx,sy,0)
pml_layers = [mp.PML(thickness=dpml,direction=mp.X)]

geometry = [ mp.Block(material=glass, size=mp.Vector3(dpml+dsub,mp.inf,mp.inf), center=mp.Vector3(-0.5*sx+0.5*(dpml+dsub),0,0)),
             mp.Block(material=glass, size=mp.Vector3(gh,gdc*gp,mp.inf), center=mp.Vector3(-0.5*sx+dpml+dsub+0.5*gh,0,0)) ]

k_point = mp.Vector3(0,0,0)

src_pos = -0.5*sx+dpml
sources = [ mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ez, center=mp.Vector3(src_pos,0,0), size=mp.Vector3(0,sy,0)) ]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    k_point=k_point,
                    sources=sources)

# fmode is the frequency corresponding to free-space wavelength lambda=0.535 um
fmode = 2.0*math.pi/0.535
xm = 0.5*sx-dpml
mflux = sim.add_eigenmode(fmode, 0, 1, mp.FluxRegion(center=mp.Vector3(xm,0,0), size=mp.Vector3(0,sy,0)))

sim.run(mp.at_beginning(mp.output_epsilon),
        until_after_sources=mp.stop_when_fields_decayed(50, mp.Ez, mp.Vector3(xm,0), 1e-6))

##################################################
# kpoint_function to compute propagation vector
# for diffraction order #mode at frequency freq
##################################################
def kpoint_func(freq, mode):
  k0 = freq
  kx = 2.0*math.pi*(mode-1)/gp
  if kx<k0:
      kz=math.sqrt(k0*k0-kx*kx)
      return mp.Vector3(kx,0.0,kz)
  return mp.Vector3(0.0, 0.0, 0.0)

##################################################
# get eigenmode coefficients of propagating modes
##################################################
for nm in range(1,10):
  kpoint = kpoint_func(fmode, nm)
  if kpoint.z > 0.0:
     alpha = sim.get_eigenmode_coefficients(mflux, [nm])
     alpha0Plus  = alpha[2*0 + 0]  # coefficient of forward-traveling fundamental mode
     alpha0Minus = alpha[2*0 + 1]  # coefficient of backward-traveling fundamental mode
     angle=math.atan2(kpoint.x,kpoint.z) * 180/math.pi
     if (mp.am_master()):
         print("freq {} order {} (angle {}) : {},{}".format(fmode,nm,angle,abs(alpha0Plus**2),abs(alpha0Minus**2)))
