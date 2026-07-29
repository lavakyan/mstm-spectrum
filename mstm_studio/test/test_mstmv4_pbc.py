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
from mstm_studio.mstm_v4 import SPR_v4
import numpy as np
import os


_in_params = {
    'wls': np.linspace(300, 800, 31),
    'n': 2.5+0.5j,
    'n_env': 1.5,
    'D': 45,  # particle size
    'W': 150}  # box size


def test_mstmv4_periodicity():
    spr = SPR_v4(_in_params['wls'])
    spr.environment_material = _in_params['n_env']
    D = _in_params['D']
    spheres = ExplicitSpheres(2, [-50, 0.0, 0.0, D/2.,
                                   50, 0.0, 0.0, D/2.],
                              mat_filename=2*[Material(_in_params['n'])])
    spr.set_spheres(spheres)
    spr.set_incident_field(fixed=True, beta_angle=0.0,
                           alpha_angle=0.0)
    # try set too small cell
    spr.set_boundary(True, D/2., _in_params['W'])
    assert(not spr._check_spheres_in_cell())
    spr.set_boundary(True, _in_params['W'], _in_params['W'])
    assert(spr._check_spheres_in_cell())
    spr.simulate()
    # spr.plot()
    assert(np.allclose(spr.transmittance,
        [0.87511, 0.88614, 0.89495, 0.90224, 0.90842, 0.91375, 0.91843,
         0.92256, 0.92626, 0.92959, 0.93261, 0.93536, 0.93789, 0.94021,
         0.94236, 0.94435, 0.9462 , 0.94792, 0.94954, 0.95105, 0.95247,
         0.95381, 0.95507, 0.95626, 0.95739, 0.95847, 0.95949, 0.96045,
         0.96138, 0.96226, 0.96309],
        rtol=1e-5))
    assert(np.allclose(spr.absorptance,
        [0.11654 , 0.10628 , 0.098103, 0.091356, 0.085653, 0.08074 ,
         0.076445, 0.07265 , 0.06926 , 0.06621 , 0.063444, 0.060925,
         0.058608, 0.056484, 0.05452 , 0.052698, 0.051002, 0.049418,
         0.047937, 0.046547, 0.04524 , 0.044009, 0.042848, 0.041746,
         0.040702, 0.03971 , 0.038769, 0.037874, 0.037016, 0.036202,
         0.035423],
        rtol=1e-5))
    assert(np.allclose(spr.reflectance,
        [0.0083471, 0.0075797, 0.0069461, 0.0064032, 0.0059281, 0.0055066,
         0.0051293, 0.0047896, 0.0044819, 0.0042023, 0.0039471, 0.0037142,
         0.0035   , 0.0033039, 0.0031232, 0.0029564, 0.0028023, 0.0026595,
         0.0025273, 0.0024043, 0.00229  , 0.0021832, 0.0020839, 0.0019908,
         0.0019038, 0.0018221, 0.0017456, 0.0016739, 0.0016059, 0.0015424,
         0.0014824],
        rtol=1e-5))
    assert(np.allclose(spr.transmittance_par,
        [0.85326, 0.86742, 0.87847, 0.88745, 0.89495, 0.90136, 0.90693,
         0.91182, 0.91617, 0.92006, 0.92358, 0.92678, 0.92971, 0.93239,
         0.93486, 0.93714, 0.93927, 0.94124, 0.94309, 0.94481, 0.94644,
         0.94796, 0.9494 , 0.95076, 0.95204, 0.95326, 0.95441, 0.95551,
         0.95656, 0.95755, 0.9585],
        rtol=1e-5))
    assert(np.allclose(spr.absorptance_par,
        [0.13657 , 0.12348 , 0.11329 , 0.10501 , 0.098114, 0.09223 ,
         0.087129, 0.082651, 0.078674, 0.075111, 0.071894, 0.068975,
         0.066293, 0.063845, 0.061586, 0.059495, 0.057552, 0.05574 ,
         0.054048, 0.052462, 0.050973, 0.049566, 0.048246, 0.046994,
         0.045808, 0.044683, 0.043617, 0.042603, 0.041631, 0.040709,
         0.039828],
        rtol=1e-5))
    assert(np.allclose(spr.reflectance_par,
        [0.010167 , 0.009097 , 0.0082467, 0.0075394, 0.0069346, 0.0064078,
         0.0059431, 0.0055295, 0.0051586, 0.0048243, 0.0045213, 0.0042463,
         0.0039946, 0.0037651, 0.0035545, 0.0033608, 0.0031823, 0.0030174,
         0.0028649, 0.0027234, 0.0025921, 0.0024695, 0.0023559, 0.0022495,
         0.0021501, 0.0020569, 0.0019698, 0.0018881, 0.0018108, 0.0017386,
         0.0016705],
        rtol=1e-5))
    assert(np.allclose(spr.transmittance_ort,
        [0.89697, 0.90486, 0.91143, 0.91704, 0.92189, 0.92615, 0.92992,
         0.9333 , 0.93635, 0.93911, 0.94163, 0.94394, 0.94607, 0.94803,
         0.94985, 0.95155, 0.95313, 0.9546 , 0.95598, 0.95728, 0.95851,
         0.95965, 0.96074, 0.96177, 0.96275, 0.96368, 0.96456, 0.96539,
         0.9662 , 0.96696, 0.96769],
        rtol=1e-5))
    assert(np.allclose(spr.absorptance_ort,
        [0.096507, 0.089077, 0.082921, 0.077697, 0.073191, 0.069249,
         0.065762, 0.062649, 0.059847, 0.057308, 0.054993, 0.052876,
         0.050922, 0.049123, 0.047454, 0.045901, 0.044452, 0.043096,
         0.041826, 0.040632, 0.039507, 0.038451, 0.03745 , 0.036498,
         0.035595, 0.034736, 0.033921, 0.033146, 0.032401, 0.031694,
         0.031017],
        rtol=1e-5))
    assert(np.allclose(spr.reflectance_ort,
        [0.0065273, 0.0060624, 0.0056455, 0.0052669, 0.0049215, 0.0046054,
         0.0043156, 0.0040496, 0.0038051, 0.0035802, 0.0033729, 0.003182 ,
         0.0030054, 0.0028426, 0.0026918, 0.002552 , 0.0024223, 0.0023017,
         0.0021897, 0.0020852, 0.0019878, 0.0018968, 0.0018119, 0.0017322,
         0.0016575, 0.0015873, 0.0015215, 0.0014597, 0.001401 , 0.0013462,
         0.0012943],
        rtol=1e-5))


def test_mstmv4_halfspace():
    spr = SPR_v4(_in_params['wls'])  # temp_dir='./temp/')
    spr.environment_material = _in_params['n_env']
    D = _in_params['D']
    spheres = ExplicitSpheres(2, [-1.05*D, 0.0, 0.0, D/2.,
                                   1.05*D, 0.0, 0.0, D/2.],
                              mat_filename=2*[Material(_in_params['n'])])
    spr.set_spheres(spheres)
    spr.set_incident_field(fixed=True)
    spr.set_boundary(False)
    spr.set_layers([Material(1.0)])
    assert(not spr._check_spheres_in_layer())
    # move spheres deeper
    spheres = ExplicitSpheres(2, [-1.05*D, 0.0, -D/1.95, D/2.,
                                   1.05*D, 0.0, -D/1.95, D/2.],
                              mat_filename=2*[Material(_in_params['n'])])
    spr.set_spheres(spheres)
    assert(spr._check_spheres_in_layer())

    spr.simulate()

    assert(np.allclose(spr.extinction,
        [1.2251 , 1.1548 , 1.0918 , 1.0348 , 0.98296, 0.93559, 0.89218,
         0.85234, 0.81566, 0.78184, 0.75056, 0.72164, 0.69474, 0.66981,
         0.64657, 0.62488, 0.6046 , 0.58559, 0.56777, 0.55101, 0.53523,
         0.52031, 0.50628, 0.49296, 0.48034, 0.46835, 0.457  , 0.4462 ,
         0.43584, 0.42603, 0.41665],
        rtol=1e-5))
    assert(np.allclose(spr.extinction_up,
        [1.422  , 1.3372 , 1.2597 , 1.1886 , 1.1234 , 1.0635 , 1.0084 ,
         0.95778, 0.91112, 0.8681 , 0.82838, 0.79172, 0.7577 , 0.72626,
         0.69704, 0.66987, 0.64455, 0.62091, 0.59885, 0.57819, 0.55882,
         0.54058, 0.52351, 0.50738, 0.49217, 0.47779, 0.46422, 0.45138,
         0.43913, 0.42757, 0.41657],
        rtol=1e-5))
    assert(np.allclose(spr.extinction_dn,
        [-1.9695e-01, -1.8240e-01, -1.6789e-01, -1.5383e-01, -1.4048e-01,
         -1.2794e-01, -1.1626e-01, -1.0544e-01, -9.5457e-02, -8.6269e-02,
         -7.7821e-02, -7.0078e-02, -6.2961e-02, -5.6448e-02, -5.0470e-02,
         -4.4985e-02, -3.9950e-02, -3.5325e-02, -3.1080e-02, -2.7175e-02,
         -2.3584e-02, -2.0269e-02, -1.7230e-02, -1.4419e-02, -1.1827e-02,
         -9.4316e-03, -7.2245e-03, -5.1859e-03, -3.2875e-03, -1.5429e-03,
          7.4897e-05],
        rtol=1e-5))
    assert(np.allclose(spr.extinction_dn_par,
        [-0.18834   , -0.17529   , -0.16195   , -0.14881   , -0.13617   ,
         -0.1242    , -0.11297   , -0.10252   , -0.092824  , -0.083874  ,
         -0.07562   , -0.068037  , -0.061053  , -0.054652  , -0.048768  ,
         -0.043363  , -0.038399  , -0.033833  , -0.029641  , -0.025783  ,
         -0.022234  , -0.018956  , -0.015951  , -0.013171  , -0.010607  ,
         -0.0082379 , -0.006055  , -0.004039  , -0.0021619 , -0.00043726,
          0.0011617 ],
        rtol=1e-5))
    assert(np.allclose(spr.scattering,
        [0.25436 , 0.22516 , 0.19898 , 0.17565 , 0.15508 , 0.13688 ,
         0.12091 , 0.10692 , 0.094688, 0.083993, 0.074637, 0.066464,
         0.059292, 0.053021, 0.04751 , 0.042662, 0.038391, 0.034617,
         0.031283, 0.028326, 0.025701, 0.023357, 0.021278, 0.019413,
         0.017744, 0.016246, 0.014902, 0.013693, 0.012594, 0.011608,
         0.010714],
        rtol=1e-5))
    assert(np.allclose(spr.scattering_up,
        [0.063372 , 0.053215 , 0.044956 , 0.038184 , 0.032597 , 0.027961 ,
         0.024092 , 0.020849 , 0.018116 , 0.015804 , 0.013838 , 0.012163 ,
         0.010725 , 0.009492 , 0.0084269, 0.0075044, 0.0067026, 0.006003 ,
         0.0053919, 0.0048552, 0.0043831, 0.0039651, 0.0035971, 0.0032694,
         0.0029779, 0.0027176, 0.0024855, 0.0022776, 0.0020896, 0.0019216,
         0.0017699],
        rtol=1e-5))
    assert(np.allclose(spr.scattering_dn,
        [0.19099  , 0.17195  , 0.15403  , 0.13747  , 0.12241  , 0.10886  ,
         0.096758 , 0.086022 , 0.076522 , 0.068142 , 0.060756 , 0.054261 ,
         0.048531 , 0.043496 , 0.039052 , 0.03513  , 0.031662 , 0.02859  ,
         0.025869 , 0.023451 , 0.021299 , 0.019375 , 0.017664 , 0.016129 ,
         0.014752 , 0.013515 , 0.012404 , 0.011403 , 0.010493 , 0.0096761,
         0.0089344],
        rtol=1e-5))
    assert(np.allclose(spr.scattering_dn_par,
        [0.2445   , 0.21553  , 0.18941  , 0.16619  , 0.14575  , 0.12787  ,
         0.11229  , 0.098763 , 0.087017 , 0.076827 , 0.067974 , 0.06029  ,
         0.053589 , 0.047761 , 0.042665 , 0.038203 , 0.034289 , 0.030845 ,
         0.027813 , 0.025132 , 0.02276  , 0.020648 , 0.01878  , 0.017108 ,
         0.015615 , 0.014278 , 0.013081 , 0.012006 , 0.01103  , 0.010157 ,
         0.0093655],
        rtol=1e-5))

    # ~ import matplotlib.pyplot as plt
    # ~ plt.plot(spr.wavelengths, spr.extinction_up, label='ext_up')
    # ~ plt.plot(spr.wavelengths, spr.extinction_dn, label='ext_dn')
    # ~ plt.plot(spr.wavelengths, spr.extinction_up_par, label='ext_up_par')
    # ~ plt.plot(spr.wavelengths, spr.extinction_dn_par, label='ext_dn_par')
    # ~ plt.legend()
    # ~ plt.show()
    # ~ plt.plot(spr.wavelengths, spr.scattering_up, label='sca_up')
    # ~ plt.plot(spr.wavelengths, spr.scattering_dn, label='sca_dn')
    # ~ plt.plot(spr.wavelengths, spr.scattering_up_par, label='sca_up_par')
    # ~ plt.plot(spr.wavelengths, spr.scattering_dn_par, label='sca_dn_par')
    # ~ plt.legend()
    # ~ plt.show()


def test_mstmv4_layer():
    ''' air pores in material slab '''
    spr = SPR_v4(_in_params['wls'], temp_dir='./temp/')
    spr.environment_material = _in_params['n_env']
    D = _in_params['D']
    spheres = ExplicitSpheres(2, [-1.05*D, 0.0, D/2. + D/10., D/2.,
                                   1.05*D, 0.0, D/2. + D/10., D/2.],
                              mat_filename=2*[Material(1.0)])
    spr.set_spheres(spheres)
    spr.set_incident_field(fixed=True)
    spr.set_boundary(False)
    spr.set_layers([Material(_in_params['n']), Material(1.0)],
                   [D/2.])
    assert(not spr._check_spheres_in_layer())
    # move spheres deeper
    spr.set_layers([Material(_in_params['n']), Material(1.0)],
                   [D + D/5.])
    assert(spr._check_spheres_in_layer())

    spr.simulate()


if __name__ == '__main__':
    # ~ test_mstmv4_periodicity()
    # ~ test_mstmv4_halfspace()
    test_mstmv4_layer()

