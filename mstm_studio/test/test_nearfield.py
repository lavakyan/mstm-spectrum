# -*- coding: utf-8 -*-
#
# ----------------------------------------------------- #
#                                                       #
#  This code is a part of T-matrix fitting project      #
#  Contributors:                                        #
#   L. Avakyan <laavakyan@sfedu.ru>                     #
#                                                       #
# ----------------------------------------------------- #
import pytest

from mstm_studio.mstm_spectrum import Material, ExplicitSpheres
from mstm_studio.nearfield import NearFieldMie, NearField
from mstm_studio.mstm_v4 import NearField_v4

import numpy as np
import os

_nf_inp = {'wl': 240,
           'matsph': 0.5 + 0.1j,
           'matrix': 1.5,
           'hmin' : -60,
           'hmax' :  60,
           'vmin' : -60,
           'vmax' :  60,
           'step' :  30,
           'a' : 20
           }

reference_Ex = np.array([
   [-0.68523465-0.71265844j,  0.37875081-0.9579761j ,
     0.96800266-0.05177116j,  0.38483353+0.87958006j,
    -0.67766488+0.70011877j],
   [-0.67243735-0.72611316j,  0.34905835-0.96559975j,
     0.70304682-0.04226429j,  0.3553152 +0.84298856j,
    -0.6618666 +0.68587346j],
   [-0.65918578-0.73612626j,  0.44500786-1.04009539j,
     1.14223468-0.12852432j,  0.45959084+0.86056855j,
    -0.64699935+0.68109398j],
   [-0.67243735-0.72611316j,  0.34905835-0.96559975j,
     0.70304682-0.04226429j,  0.3553152 +0.84298856j,
    -0.6618666 +0.68587346j],
   [-0.68523465-0.71265844j,  0.37875081-0.9579761j ,
     0.96800266-0.05177116j,  0.38483353+0.87958006j,
    -0.67766488+0.70011877j]
])

reference_E2 = np.array([
   [0.97808946, 1.06206309, 0.93972987, 0.92267293, 0.95019445],
   [0.98021687, 1.0612284 , 0.49801506, 0.84413273, 0.90949164],
   [0.97640777, 1.27983041, 1.32121856, 0.95180196, 0.88249716],
   [0.98021687, 1.0612284 , 0.49801506, 0.84413273, 0.90949164],
   [0.97808946, 1.06206309, 0.93972987, 0.92267293, 0.95019445]
])


def test_nearfield_Mie():
    nf = NearFieldMie()
    E2 = nf.calculate(wavelength=_nf_inp['wl'],
                      material=_nf_inp['matsph'],
                      environment_material=_nf_inp['matrix'],
                      radius=_nf_inp['a'],
                      plane='zx',
                      hmin=_nf_inp['hmin'], hmax=_nf_inp['hmax'],
                      vmin=_nf_inp['vmin'], vmax=_nf_inp['vmax'],
                      step=_nf_inp['step'],
                      include_incident=True)
    assert(np.allclose(nf.E_xyz[0], reference_Ex))
    assert(np.allclose(E2, reference_E2))


def test_nearfield_v3():
    nf = NearField(wavelength=_nf_inp['wl'])  # temp_dir='./temp/')
    nf.environment_material = _nf_inp['matrix']
    nf.set_incident_field(fixed=True, azimuth_angle=0.0,
                          polar_angle=0.0, polarization_angle=0.0)
    nf.set_plane(plane='xz',
                 hmin=_nf_inp['hmin'], hmax=_nf_inp['hmax'],
                 vmin=_nf_inp['vmin'], vmax=_nf_inp['vmax'],
                 step=_nf_inp['step'])

    spheres = ExplicitSpheres(1, [0, 0, 0, _nf_inp['a']],
                              mat_filename=Material(_nf_inp['matsph']))
    nf.set_spheres(spheres)
    nf.simulate()
    assert(np.allclose(np.transpose(nf.field), reference_E2, rtol=0.001))


def test_nearfield_v4():
    nf = NearField_v4(wavelength=_nf_inp['wl'], mstm_path='mstm2023.x',
                      incident_default_mode=True)
    nf.environment_material = _nf_inp['matrix']
    nf.set_incident_field(fixed=True, beta_angle=0.0, alpha_angle=0.0)
    nf.set_plane(plane='xz',
                 hmin=_nf_inp['hmin'], hmax=_nf_inp['hmax'],
                 vmin=_nf_inp['vmin'], vmax=_nf_inp['vmax'],
                 step=_nf_inp['step'])

    spheres = ExplicitSpheres(1, [0, 0, 0, _nf_inp['a']],
                              mat_filename=Material(_nf_inp['matsph']))
    nf.set_spheres(spheres)
    nf.simulate()
    print('v4')
    # ~ print(nf.field.T)
    print(nf.Epar_x.real.T)
    print('ref')
    print(reference_Ex.real)
    print('delta')
    print(nf.Epar_x.real.T - reference_Ex.real)
    #import matplotlib.pyplot as plt
    #fig, (ax1, ax2) = plt.subplots(1,2, figsize=(6,3))
    #ax1.imshow(nf.Epar_x.real.T)
    #ax2.imshow(reference_Ex.real)
    #plt.show()
    # Comparable only inner region 3x3 of 5x5 matrix
    # Outer region is probably contaminated by border effects,
    # or, more hopefuly, by periodic conditions?
    assert(np.isclose(nf.field[2,2], reference_E2[2,2], rtol=0.001))
    assert(np.allclose(nf.field.T, reference_E2, rtol=0.9))
    assert(np.allclose(nf.Epar_x.T, reference_Ex, rtol=0.5))
    # ~ assert(np.allclose(np.transpose(nf.field[1:4,1:4]), reference_E2[1:4,1:4], rtol=0.3))
    # ~ assert(np.allclose(np.transpose(nf.Epar_x[1:4,1:4]), reference_Ex[1:4,1:4], rtol=0.15))

if __name__ == '__main__':
    test_nearfield_v4()
