#! /usr/bin/env python
import matplotlib.pyplot as plt
from mstm_studio.mstm_v4 import NearField_v4
from mstm_studio.mstm_spectrum import Material, ExplicitSpheres

wl = 300  # wavelength, nm
matsph = 1.34 + 0.964j  # silver at 300 nm
a = 5  # particle radii, nm
layers_n = [1.5, 1.33, 1.0]  # refr. indeces of layers
depth = 20  # depth of 2d layer
# near field calculation region:
hmin, hmax, vmin, vmax, step = -15, 15, -10, 40, 0.1


nf = NearField_v4(wavelength=wl)
nf.environment_material = Material(layers_n[0])
spheres = ExplicitSpheres(1, [0, 0, a+5, a],
                          mat_filename=Material(matsph))
nf.set_plane(plane='xz', hmin=hmin, hmax=hmax,
             vmin=vmin, vmax=vmax, step=step, offset=0)
nf.set_layers([Material(layers_n[1]),
               Material(layers_n[2])],
              [depth])
nf.set_spheres(spheres)

nf.simulate()  # do calc

# do plot
fig, ax = plt.subplots(1, 1, figsize=(6, 6))
nf.plot(fig=fig, axs=ax)
ax.text(0, a+5, 'Ag', color='white')
ax.text(-14, -5, f'n={layers_n[0]:.2f}', color='white')
ax.text(-14, 0.5*depth, f'n={layers_n[1]:.2f}', color='white')
ax.text(-14, 1.5*depth, f'n={layers_n[2]:.2f}', color='white')
plt.savefig('example_silver_in_layer.png')
plt.show()
