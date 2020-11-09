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

from mstm_studio import mstm_spectrum
import numpy as np
import os

# ~ import numpy as np
# ~ from numpy.random import lognormal
# ~ from scipy import interpolate
# ~ import subprocess
# ~ import os   # to delete files after calc.
# ~ import sys  # to check whether running on Linux or Windows
# ~ import datetime
# ~ import time
# ~ import tempfile  # to run mstm in temporary directory
# ~ try:
    # ~ import matplotlib.pyplot as plt
# ~ except ImportError:
    # ~ pass
# ~ # use input in both python2 and python3
# ~ try:
    # ~ input = raw_input
# ~ except NameError:
    # ~ pass
# ~ # use xrange in both python2 and python3
# ~ try:
    # ~ xrange
# ~ except NameError:
    # ~ xrange = range

def test_spheres():
    print('Test Spheres')
    print('Overlap tests')
    spheres = mstm_spectrum.Spheres()
    print('  Test not overlapped... ')
    spheres.x = [-5, 5]
    spheres.y = [0, 0]
    spheres.z = [0, 0]
    spheres.a = [4, 4]
    assert(not spheres.check_overlap())
    print('  Test overlapped... ')
    spheres.a = [5, 5]
    assert(spheres.check_overlap())
    print('  Test nested... ')
    spheres.x = [0, 0]
    spheres.a = [2, 5]
    assert(not spheres.check_overlap())
    spheres.a = [5, 3]
    assert(not spheres.check_overlap())


def test_materials():
    print('Test Materials')
    mat = mstm_spectrum.Material(os.path.join('..', 'nk', 'etaGold.txt'))
    # mat.plot()
    mat1 = mstm_spectrum.Material(os.path.join('..', 'nk', 'etaSilver.txt'))
    mat3 = mstm_spectrum.Material('glass')
    mat5 = mstm_spectrum.Material(1.5)
    mat6 = mstm_spectrum.Material('2.0+0.5j')
    mat7 = mstm_spectrum.Material('mat7', wls=np.linspace(300, 800, 100),
                    nk=np.linspace(-10, 5, 100) + 1j * np.linspace(0, 10, 100))
    mat8 = mstm_spectrum.Material('mat7', wls=np.linspace(300, 800, 100),
                    eps=np.linspace(-10, 5, 100) + 1j * np.linspace(0, 10, 100))
    assert(np.isclose(mat.get_n(800), 0.15436829401))  # etaGold
    assert(np.isclose(mat1.get_n(800), 0.03604950826))  # etaSilver
    assert(np.isclose(mat3.get_n(800), 1.66))  # Glass (constant)
    assert(np.isclose(mat3.get_k(800), 0.00))
    assert(np.isclose(mat5.get_n(550), 1.5))  # n=1.5 material
    assert(np.isclose(mat6.get_n(550), 2.0))  # n=2.0+0.5j material
    assert(np.isclose(mat6.get_k(550), 0.5))
    assert(np.isclose(mat7.get_n(550), -2.5))  # nk material
    assert(np.isclose(mat7.get_k(550), 5.0))
    assert(np.isclose(mat8.get_n(550), 1.243014470))  # eps material
    assert(np.isclose(mat8.get_k(550), 2.011239667))


def test_SPR():
    wls = np.linspace(300, 800, 100)
    # create SPR object
    spr = mstm_spectrum.SPR(wls)
    spr.environment_material = 'glass'
    # spr.set_spheres(SingleSphere(0.0, 0.0, 0.0, 25.0, 'etaGold.txt'))
    spheres = mstm_spectrum.ExplicitSpheres(2, [0, 0, 0, 10, 0, 0, 0, 12],
                              mat_filename=['../nk/etaGold.txt',
                                            '../nk/etaSilver.txt'])
    # spheres = ExplicitSpheres(2, [0,0,0,20,0,0,0,21],
    #                           mat_filename='etaGold.txt')
    spr.set_spheres(spheres)
    # spr.set_spheres(LogNormalSpheres(27, 0.020, 0.9, 0.050 ))
    # calculate!
    # spr.command = ''
    _, exts = spr.simulate()
    assert(np.allclose(
        (exts[20], exts[50], exts[-20]), (3.3279, 1.9906, 0.056079)))

