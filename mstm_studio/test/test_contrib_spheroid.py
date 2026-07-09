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

from mstm_studio.mstm_spectrum import Material
from mstm_studio.contrib_spheroid import SpheroidSP
from mstm_studio.contributions import MieSingleSphere

import numpy as np
import os

def test_spheroidSP():
    n_env = 1.5
    mat_gold = Material(os.path.join('..', 'nk', 'etaGold.txt'))
    wls = np.linspace(300, 800, 15)
    npsize = 100  # diameter of nanoparticle
    sph = SpheroidSP(wavelengths=wls)
    sph.set_material(mat_gold, n_env)
    ext_sph = sph.calculate([1, npsize, 1.0])

    mie = MieSingleSphere(name='mie', wavelengths=wls)
    mie.set_material(mat_gold, n_env)
    ext_mie = mie.calculate([1, npsize])

    # import matplotlib.pyplot as plt
    # plt.plot(wls, ext_sph, label='sph')
    # plt.plot(wls, ext_mie, label='mie')
    # plt.legend()
    # plt.show()
    assert(np.allclose(ext_sph, ext_mie))


