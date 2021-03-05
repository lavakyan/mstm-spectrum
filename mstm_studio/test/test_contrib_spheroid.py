# -*- coding: utf-8 -*-
#
# ----------------------------------------------------- #
#                                                       #
#  This code is a part of T-matrix fitting project      #
#  Contributors:                                        #
#   L. Avakyan <laavakyan@sfedu.ru>                     #
#                                                       #
# ----------------------------------------------------- #
from __future__ import print_function
from __future__ import division

import pytest

from mstm_studio.mstm_spectrum import Material
from mstm_studio.contrib_spheroid import SpheroidSP
from mstm_studio.contributions import MieSingleSphere

import numpy as np
import os


def test_spheroidSP():
    mat_gold = Material(os.path.join('..', 'nk', 'etaGold.txt'))
    wls = np.linspace(300, 800, 51)
    npsize = 10  # diameter of nanoparticle
    sph = SpheroidSP(wavelengths=wls)
    sph.set_material(mat_gold, 1.5)
    ext_sph = sph.calculate([1, npsize, 1.0])

    mie = MieSingleSphere(name='mie', wavelengths=wls)
    mie.set_material(mat_gold, 1.5)
    ext_mie = mie.calculate([1, npsize])

    assert(np.allclose(ext_sph, ext_mie))

