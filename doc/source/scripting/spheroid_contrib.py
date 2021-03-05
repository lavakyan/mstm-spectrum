
from mstm_studio.alloy_AuAg import AlloyAuAg
from mstm_studio.contrib_spheroid import SpheroidSP
import numpy as np


wls = np.linspace(300, 800, 51)  # range for calculation, in nm
SIZE = 20  # nm, particle diameter
ASPECT = 1.5  # a / c = horiz. axis / rot. axis

sph = SpheroidSP(wavelengths=wls)  # create object
sph.set_material(AlloyAuAg(x_Au=1), 1.5)  # particle and matrix refr. ind.

fig, axs = sph.plot_shape([1, SIZE, ASPECT])
fig.savefig('spheroid_shape.png', bbox_inches='tight')

ext_sph = sph.calculate([1, SIZE, ASPECT])
fig, axs = sph.plot([1, SIZE, ASPECT])  # scale, diameter, aspect
fig.savefig('spheroid_ext.png', bbox_inches='tight')
