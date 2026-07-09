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
import numpy as np
from mstm_studio.contributions import MieSingleSphere
from mstm_studio.mstm_spectrum import SPR, Material, SingleSphere
from mstm_studio.mstm_v4 import SPR_v4


_spec_in = {  # input data for the tests
    'wls' : np.linspace(300, 700, 10),
    'a' : 100,
    'n' : 2.5+0.5j,
    'nenv' : 1.33,
    }

reference_ext = np.array(
    [4.09671076, 4.21590039, 4.3395197 , 4.18193652, 3.97191584,
     3.83814765, 3.68796047, 3.45511005, 3.15756782, 2.84275593])
reference_sca = np.array(
    [2.0177035 , 2.17920573, 2.25789734, 2.18030795, 2.08751548,
     1.9888349 , 1.84504133, 1.65904351, 1.45945261, 1.26798841])
reference_abs = np.array(
    [2.07900725, 2.03669465, 2.08162236, 2.00162857, 1.88440036,
     1.84931275, 1.84291914, 1.79606654, 1.69811522, 1.57476752])


def test_spec_Mie():
    mie = MieSingleSphere(name='mie', wavelengths=_spec_in['wls'])
    mie.set_material(_spec_in['n'], _spec_in['nenv'])
    mie.calculate([1, 2*_spec_in['a']])

    assert(np.allclose(mie.qext, reference_ext))
    assert(np.allclose(mie.qsca, reference_sca))
    assert(np.allclose(mie.qabs, reference_abs))


def test_spec_SPRv3():
    spr = SPR(_spec_in['wls'])
    spr.environment_material = _spec_in['nenv']
    spr.set_spheres(SingleSphere(0.0, 0.0, 0.0,
                                 _spec_in['a'],
                                 Material(_spec_in['n'])))
    spr.set_incident_field(fixed=False)
    spr.simulate()

    assert(np.allclose(spr.extinction, reference_ext, rtol=1e-4))
    assert(np.allclose(spr.scattering, reference_sca, rtol=1e-3))
    assert(np.allclose(spr.absorbtion, reference_abs, rtol=1e-4))


def test_spec_SPRv4():
    spr = SPR_v4(_spec_in['wls'], mstm_path='~/bin/mstm2023.x')
                 #temp_dir='./temp/')
    spr.environment_material = _spec_in['nenv']
    spr.set_spheres(SingleSphere(0.0, 0.0, 0.0,
                                 _spec_in['a'],
                                 Material(_spec_in['n'])))
    spr.set_incident_field(fixed=False)
    spr.simulate()

    assert(np.allclose(spr.extinction, reference_ext, rtol=1e-4))
    assert(np.allclose(spr.scattering, reference_sca, rtol=1e-3))
    assert(np.allclose(spr.absorbtion, reference_abs, rtol=1e-4))
