
from mstm_studio.nearfield import NearField
from mstm_studio.mstm_spectrum import ExplicitSpheres
from mstm_studio.alloy_AuAg import AlloyAuAg

mat = AlloyAuAg(x_Au=0.0)       # silver material
nf = NearField(wavelength=385)  # near the resonance
nf.environment_material = 1.5   # glass matrix
nf.set_plane(plane='xz', hmin=-10, hmax=20, vmin=-15, vmax=15, step=0.25)

spheres = ExplicitSpheres(2, [0, 0, 0, 5, 0, 0, 11, 3],
                          mat_filename=2*[mat])
nf.set_spheres(spheres)

nf.simulate()

nf.plot()
