from mstm_studio.mstm_spectrum import LogNormalSpheres
from mstm_studio.alloy_AuAg import AlloyAuAg

spheres = LogNormalSpheres(9, 10.0, 0.2, 5., AlloyAuAg(1.))
while spheres.check_overlap():
    print('spheres are overlapping, regenerating...')
    spheres = LogNormalSpheres(9, 10.0, 0.2, 5., AlloyAuAg(1.))
print(spheres.a)


# ~ Box size estimated as: 77.0 nm
# ~ Desired number of particles: 9
# ~ Number of particles in a box: 8
# ~ Resulted number of particles: 8
# ~ spheres are overlapping, regenerating...
# ~ Box size estimated as: 77.0 nm
# ~ Desired number of particles: 9
# ~ Number of particles in a box: 8
# ~ Resulted number of particles: 8
# ~ [ 9.307773    8.61185299  9.92867988  8.84140858  9.87175352  8.71090184
  # ~ 9.71505038 12.40459688]
