#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from mstm_studio.mstm_v4 import SPR_v4
from mstm_studio.mstm_spectrum import Material, Spheres, SingleSphere

# specify the path to gold dielectric function
mat1 = Material('nk/etaGold.txt')

# calculation parameters
wls = np.linspace(400, 700, 61)  # wavelength, nm
cell = 30  # distance between particle centers, nm
a = 10  # particle radius, nm

# create caclulator object
spr = SPR_v4(wls, mstm_path='~/bin/mstm2023.x')  # , temp_dir='./temp/')
spr.environment_material = 1.5  # glass

# 1. isolated particle
spheres = SingleSphere(0, 0, 0, a, mat1)
spr.set_spheres(spheres)
spr.set_incident_field(fixed=True, beta_angle=0.0, alpha_angle=0.0)
spr.simulate()
isol_ext = spr.extinction_par

# 2. N spheres chain
def ext_of_spheres_chain(N):
    ''' helper function to calc
    extinction for chain '''
    spheres = Spheres()  # empty
    for i in range(-N//2, N//2+1):
        spheres.append(
            SingleSphere(x=i*cell, y=0, z=0, a=a,
                         mat_filename=mat1))
    spr.set_spheres(spheres)
    spr.simulate()
    return spr.extinction_par

chain3_ext = ext_of_spheres_chain(N=3)
chain5_ext = ext_of_spheres_chain(N=5)
chain7_ext = ext_of_spheres_chain(N=7)

# use of PBC
spheres = SingleSphere(0, 0, 0, a, mat1)
spr.set_spheres(spheres)
spr.set_boundary(pbc=True,  # use periodic boundaries
                 cell_x=cell,  # chain direction
                 cell_y=100)  # some gap in off-chain direction
                              # tests show that effic. are indep. of it
spr.simulate()
pbc_ext = spr.extinction_par

# do plots:
plt.figure(figsize=(4,3))
plt.plot(wls, chain7_ext, linewidth=1, linestyle=':', color='green', label='chain of 7')
plt.plot(wls, chain5_ext, 'g-.', linewidth=1, label='chain of 5')
plt.plot(wls, chain3_ext, 'g-',  linewidth=1, label='chain of 3')
plt.plot(wls, pbc_ext, 'b-', linewidth=1.5, label='PBC')
plt.plot(wls, isol_ext, 'r--', linewidth=1.5, label='isolated')
plt.xlabel(r'Wavelength, nm')
plt.ylabel(r'Ext. effic.')
plt.legend()
plt.tight_layout()
plt.savefig('example_periodic_chain.png')
plt.show()
