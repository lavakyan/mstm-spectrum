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

from mstm_studio.mstm_spectrum import Material, SPR,\
                                      ExplicitSpheres, Spheres
from mstm_studio.mstm_v4 import SPR_v4
import numpy as np
import os


def test_spheres():
    print('Test Spheres')
    print('Overlap tests')
    spheres = Spheres()
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
    mat = Material(os.path.join('..', 'nk', 'etaGold.txt'))
    # mat.plot()
    mat1 = Material(os.path.join('..', 'nk', 'etaSilver.txt'))
    mat3 = Material('glass')
    mat5 = Material(1.5)
    mat6 = Material('2.0+0.5j')
    mat7 = Material('mat7', wls=np.linspace(300, 800, 100),
                    nk=np.linspace(-10, 5, 100) + 1j * np.linspace(0, 10, 100))
    mat8 = Material('mat7', wls=np.linspace(300, 800, 100),
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


def test_mstmv3_AuAg_coreshell():
    ref_exts = \
      [1.7283  , 1.5937  , 1.4142  , 1.1876  , 1.0103  , 1.0952  ,
       1.6424  , 2.468   , 2.7153  , 2.8001  , 2.8943  , 2.9381  ,
       2.9432  , 2.941   , 2.9477  , 2.9556  , 2.9537  , 2.9325  ,
       2.895   , 2.8566  , 2.8272  , 2.8054  , 2.7869  , 2.7644  ,
       2.7397  , 2.7187  , 2.709   , 2.7135  , 2.7288  , 2.7446  ,
       2.7605  , 2.7671  , 2.7767  , 2.7945  , 2.835   , 2.9036  ,
       3.0047  , 3.1338  , 3.2812  , 3.4285  , 3.5205  , 3.4874  ,
       3.2755  , 2.8963  , 2.4363  , 1.9922  , 1.6195  , 1.3231  ,
       1.0927  , 0.91192 , 0.76418 , 0.64392 , 0.54536 , 0.4609  ,
       0.39244 , 0.3356  , 0.28919 , 0.25042 , 0.21966 , 0.19321 ,
       0.17134 , 0.15279 , 0.13723 , 0.12292 , 0.11049 , 0.099222,
       0.089345, 0.0806  , 0.072448, 0.065765, 0.059703, 0.054431,
       0.050077, 0.046343, 0.043312, 0.040489, 0.038094, 0.036043,
       0.034298, 0.032521, 0.031028, 0.029624, 0.028304, 0.026845,
       0.025694, 0.02461 , 0.023595, 0.022632, 0.021738, 0.02088 ,
       0.020261, 0.019571, 0.019004, 0.018444, 0.01793 , 0.017496,
       0.017003, 0.0166  , 0.016215, 0.015664]

    wls = np.linspace(300, 800, 100)
    spr = SPR(wls)
    spr.environment_material = 1.5
    spheres = ExplicitSpheres(2, [0, 0, 0, 10, 0, 0, 0, 12],
        mat_filename=[os.path.join('..', 'nk', 'etaGold.txt'),
                      os.path.join('..', 'nk', 'etaSilver.txt')])
    spr.set_spheres(spheres)
    _, exts = spr.simulate()
    assert(np.allclose(exts, ref_exts))


_in_2sph = {
    'wls': np.linspace(300, 800, 51),
    'n': 2.5+0.5j,
    'n_env': 1.5,
    'D': 100}

ref_ext_par = np.array([
       4.3226 , 4.1866 , 4.0495 , 3.9143 , 3.7819 , 3.652  , 3.5236 ,
       3.3964 , 3.2695 , 3.1435 , 3.0189 , 2.8965 , 2.7771 , 2.6614 ,
       2.5501 , 2.4438 , 2.3423 , 2.2462 , 2.1557 , 2.0704 , 1.9899 ,
       1.9147 , 1.8441 , 1.7778 , 1.7159 , 1.6575 , 1.6028 , 1.5514 ,
       1.5031 , 1.4575 , 1.4147 , 1.3742 , 1.336  , 1.2998 , 1.2654 ,
       1.2329 , 1.2021 , 1.1727 , 1.1447 , 1.118  , 1.0925 , 1.0683 ,
       1.045  , 1.0226 , 1.0012 , 0.98076, 0.96092, 0.94205, 0.92378,
       0.90624, 0.88943])
ref_abs_par = np.array([
       2.3607 , 2.3106 , 2.259  , 2.2086 , 2.1606 , 2.1148 , 2.0703 ,
       2.0264 , 1.982  , 1.9368 , 1.8907 , 1.8438 , 1.7964 , 1.7487 ,
       1.7013 , 1.6543 , 1.6081 , 1.5629 , 1.5191 , 1.4766 , 1.4355 ,
       1.3962 , 1.3583 , 1.3221 , 1.2876 , 1.2545 , 1.2229 , 1.1928 ,
       1.1641 , 1.1366 , 1.1105 , 1.0855 , 1.0617 , 1.0388 , 1.017  ,
       0.9961 , 0.97615, 0.95698, 0.93856, 0.92087, 0.90389, 0.8876 ,
       0.87186, 0.85664, 0.84207, 0.828  , 0.81429, 0.80119, 0.78843,
       0.77612, 0.76426])
ref_sca_par = np.array([
       1.9619 , 1.876  , 1.7905 , 1.7057 , 1.6213 , 1.5372 , 1.4533 ,
       1.37   , 1.2875 , 1.2067 , 1.1281 , 1.0527 , 0.98068, 0.91263,
       0.84883, 0.78942, 0.73421, 0.6833 , 0.63663, 0.59377, 0.55443,
       0.51858, 0.48573, 0.45568, 0.42824, 0.40299, 0.37985, 0.35857,
       0.33899, 0.32087, 0.3042 , 0.28869, 0.27433, 0.26094, 0.24846,
       0.23682, 0.22596, 0.21576, 0.20616, 0.19714, 0.18865, 0.18067,
       0.1731 , 0.16592, 0.15916, 0.15276, 0.14663, 0.14087, 0.13535,
       0.13012, 0.12517])
ref_ext_ort = np.array([
       4.1635 , 4.0021 , 3.8564 , 3.7271 , 3.6132 , 3.5116 , 3.4193 ,
       3.333  , 3.2502 , 3.1685 , 3.0857 , 3.0008 , 2.9127 , 2.8212 ,
       2.7268 , 2.6299 , 2.5313 , 2.4323 , 2.3341 , 2.2375 , 2.1432 ,
       2.0525 , 1.9655 , 1.8827 , 1.8045 , 1.7306 , 1.6612 , 1.5962 ,
       1.5355 , 1.4787 , 1.4259 , 1.3765 , 1.3304 , 1.2873 , 1.2471 ,
       1.2094 , 1.1743 , 1.1413 , 1.1102 , 1.0809 , 1.0534 , 1.0275 ,
       1.003  , 0.97965, 0.95767, 0.93679, 0.91676, 0.89788, 0.87975,
       0.86249, 0.84607])
ref_abs_ort = np.array([
       2.305  , 2.2414 , 2.1773 , 2.116  , 2.0594 , 2.0087 , 1.9642 ,
       1.9256 , 1.8925 , 1.8637 , 1.838  , 1.8142 , 1.7909 , 1.7672 ,
       1.7422 , 1.7154 , 1.6864 , 1.6551 , 1.6219 , 1.5868 , 1.5503 ,
       1.5128 , 1.4748 , 1.4365 , 1.3986 , 1.361  , 1.3242 , 1.2883 ,
       1.2536 , 1.2199 , 1.1876 , 1.1565 , 1.1268 , 1.0982 , 1.071  ,
       1.0449 , 1.0202 , 0.99647, 0.97381, 0.95218, 0.93154, 0.91188,
       0.89301, 0.87491, 0.8577 , 0.8412 , 0.82525, 0.81011, 0.79548,
       0.78147, 0.76807])
ref_sca_ort = np.array([
       1.8585  , 1.7607  , 1.679   , 1.6111  , 1.5538  , 1.5029  ,
       1.4551  , 1.4074  , 1.3578  , 1.3048  , 1.2477  , 1.1866  ,
       1.1218  , 1.054   , 0.98452 , 0.91449 , 0.84492 , 0.77712 ,
       0.71222 , 0.65067 , 0.59293 , 0.53968 , 0.49072 , 0.44616 ,
       0.40594 , 0.36956 , 0.33702 , 0.30791 , 0.28197 , 0.25881 ,
       0.23829 , 0.21994 , 0.20364 , 0.18909 , 0.1761  , 0.16449 ,
       0.15413 , 0.14479 , 0.13637 , 0.12877 , 0.12189 , 0.11566 ,
       0.10996 , 0.10473 , 0.099973, 0.095588, 0.091505, 0.087765,
       0.084266, 0.081018, 0.077997])


def test_mstmv3_2spheres():
    spr = SPR(_in_2sph['wls'])
    spr.environment_material = _in_2sph['n_env']
    D = _in_2sph['D']
    spheres = ExplicitSpheres(2, [-1.5*D, 0.0, 0.0, D/2.,
                                   1.5*D, 0.0, 0.0, D/2.],
        mat_filename=2*[Material(_in_2sph['n'])])
    spr.set_spheres(spheres)
    spr.set_incident_field(fixed=True, azimuth_angle=0.0,
                           polar_angle=0.0, polarization_angle=0.0)
    spr.simulate()
    assert(np.allclose(spr.extinction_par, ref_ext_par))
    assert(np.allclose(spr.absorbtion_par, ref_abs_par))
    assert(np.allclose(spr.scattering_par, ref_sca_par))
    assert(np.allclose(spr.extinction_ort, ref_ext_ort))
    assert(np.allclose(spr.absorbtion_ort, ref_abs_ort))
    assert(np.allclose(spr.scattering_ort, ref_sca_ort))


def test_mstmv4_2spheres():
    spr = SPR_v4(_in_2sph['wls'], mstm_path='mstm2023.x')
    spr.environment_material = _in_2sph['n_env']
    D = _in_2sph['D']
    spheres = ExplicitSpheres(2, [-1.5*D, 0.0, 0.0, D/2.,
                                   1.5*D, 0.0, 0.0, D/2.],
        mat_filename=2*[Material(_in_2sph['n'])])
    spr.set_spheres(spheres)
    spr.set_incident_field(fixed=True, beta_angle=0.0,
                           alpha_angle=0.0)
    spr.simulate()
    assert(np.allclose(spr.extinction_par, ref_ext_par, rtol=1e-5))
    assert(np.allclose(spr.absorbtion_par, ref_abs_par, rtol=1e-4))
    assert(np.allclose(spr.scattering_par, ref_sca_par, rtol=1e-4))
    assert(np.allclose(spr.extinction_ort, ref_ext_ort, rtol=1e-4))
    assert(np.allclose(spr.absorbtion_ort, ref_abs_ort, rtol=1e-3))
    assert(np.allclose(spr.scattering_ort, ref_sca_ort, rtol=1e-3))


if __name__ == '__main__':
    test_mstmv3_2spheres()
    test_mstmv4_2spheres()
