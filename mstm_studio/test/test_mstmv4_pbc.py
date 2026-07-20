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


_in_periodic = {
    'wls': np.linspace(300, 800, 31),
    'n': 2.5+0.5j,
    'n_env': 1.5,
    'D': 45,  # particle size
    'W': 150}  # box size

ref_ext = np.array([])


ref_T = np.array(
      [0.87511, 0.88614, 0.89495, 0.90224, 0.90842, 0.91375, 0.91843,
       0.92256, 0.92626, 0.92959, 0.93261, 0.93536, 0.93789, 0.94021,
       0.94236, 0.94435, 0.9462 , 0.94792, 0.94954, 0.95105, 0.95247,
       0.95381, 0.95507, 0.95626, 0.95739, 0.95847, 0.95949, 0.96045,
       0.96138, 0.96226, 0.96309])
ref_A = np.array(
      [0.11654 , 0.10628 , 0.098103, 0.091356, 0.085653, 0.08074 ,
       0.076445, 0.07265 , 0.06926 , 0.06621 , 0.063444, 0.060925,
       0.058608, 0.056484, 0.05452 , 0.052698, 0.051002, 0.049418,
       0.047937, 0.046547, 0.04524 , 0.044009, 0.042848, 0.041746,
       0.040702, 0.03971 , 0.038769, 0.037874, 0.037016, 0.036202,
       0.035423])
ref_R = np.array(
      [0.0083471, 0.0075797, 0.0069461, 0.0064032, 0.0059281, 0.0055066,
       0.0051293, 0.0047896, 0.0044819, 0.0042023, 0.0039471, 0.0037142,
       0.0035   , 0.0033039, 0.0031232, 0.0029564, 0.0028023, 0.0026595,
       0.0025273, 0.0024043, 0.00229  , 0.0021832, 0.0020839, 0.0019908,
       0.0019038, 0.0018221, 0.0017456, 0.0016739, 0.0016059, 0.0015424,
       0.0014824])
ref_T_par = np.array([
       0.85326, 0.86742, 0.87847, 0.88745, 0.89495, 0.90136, 0.90693,
       0.91182, 0.91617, 0.92006, 0.92358, 0.92678, 0.92971, 0.93239,
       0.93486, 0.93714, 0.93927, 0.94124, 0.94309, 0.94481, 0.94644,
       0.94796, 0.9494 , 0.95076, 0.95204, 0.95326, 0.95441, 0.95551,
       0.95656, 0.95755, 0.9585 ])
ref_A_par = np.array([
       0.13657 , 0.12348 , 0.11329 , 0.10501 , 0.098114, 0.09223 ,
       0.087129, 0.082651, 0.078674, 0.075111, 0.071894, 0.068975,
       0.066293, 0.063845, 0.061586, 0.059495, 0.057552, 0.05574 ,
       0.054048, 0.052462, 0.050973, 0.049566, 0.048246, 0.046994,
       0.045808, 0.044683, 0.043617, 0.042603, 0.041631, 0.040709,
       0.039828])
ref_R_par = np.array(
      [0.010167 , 0.009097 , 0.0082467, 0.0075394, 0.0069346, 0.0064078,
       0.0059431, 0.0055295, 0.0051586, 0.0048243, 0.0045213, 0.0042463,
       0.0039946, 0.0037651, 0.0035545, 0.0033608, 0.0031823, 0.0030174,
       0.0028649, 0.0027234, 0.0025921, 0.0024695, 0.0023559, 0.0022495,
       0.0021501, 0.0020569, 0.0019698, 0.0018881, 0.0018108, 0.0017386,
       0.0016705])
ref_T_ort = np.array([
       0.89697, 0.90486, 0.91143, 0.91704, 0.92189, 0.92615, 0.92992,
       0.9333 , 0.93635, 0.93911, 0.94163, 0.94394, 0.94607, 0.94803,
       0.94985, 0.95155, 0.95313, 0.9546 , 0.95598, 0.95728, 0.95851,
       0.95965, 0.96074, 0.96177, 0.96275, 0.96368, 0.96456, 0.96539,
       0.9662 , 0.96696, 0.96769])
ref_A_ort = np.array([
       0.096507, 0.089077, 0.082921, 0.077697, 0.073191, 0.069249,
       0.065762, 0.062649, 0.059847, 0.057308, 0.054993, 0.052876,
       0.050922, 0.049123, 0.047454, 0.045901, 0.044452, 0.043096,
       0.041826, 0.040632, 0.039507, 0.038451, 0.03745 , 0.036498,
       0.035595, 0.034736, 0.033921, 0.033146, 0.032401, 0.031694,
       0.031017])
ref_R_ort = np.array(
      [0.0065273, 0.0060624, 0.0056455, 0.0052669, 0.0049215, 0.0046054,
       0.0043156, 0.0040496, 0.0038051, 0.0035802, 0.0033729, 0.003182 ,
       0.0030054, 0.0028426, 0.0026918, 0.002552 , 0.0024223, 0.0023017,
       0.0021897, 0.0020852, 0.0019878, 0.0018968, 0.0018119, 0.0017322,
       0.0016575, 0.0015873, 0.0015215, 0.0014597, 0.001401 , 0.0013462,
       0.0012943])


def test_mstmv4_periodicity():
    spr = SPR_v4(_in_periodic['wls'])
    spr.environment_material = _in_periodic['n_env']
    D = _in_periodic['D']
    spheres = ExplicitSpheres(2, [-50, 0.0, 0.0, D/2.,
                                   50, 0.0, 0.0, D/2.],
                              mat_filename=2*[Material(_in_periodic['n'])])
    spr.set_spheres(spheres)
    spr.set_incident_field(fixed=True, beta_angle=0.0,
                           alpha_angle=0.0)
    # try set too small cell
    spr.set_boundary(True, D/2., _in_periodic['W'])
    assert(not spr._check_spheres_in_cell())
    spr.set_boundary(True, _in_periodic['W'], _in_periodic['W'])
    assert(spr._check_spheres_in_cell())
    spr.simulate()
    # spr.plot()
    assert(np.allclose(spr.transmittance, ref_T, rtol=1e-5))
    assert(np.allclose(spr.absorptance, ref_A, rtol=1e-5))
    assert(np.allclose(spr.reflectance, ref_R, rtol=1e-5))
    assert(np.allclose(spr.transmittance_par, ref_T_par, rtol=1e-5))
    assert(np.allclose(spr.absorptance_par, ref_A_par, rtol=1e-5))
    assert(np.allclose(spr.reflectance_par, ref_R_par, rtol=1e-5))
    assert(np.allclose(spr.transmittance_ort, ref_T_ort, rtol=1e-5))
    assert(np.allclose(spr.absorptance_ort, ref_A_ort, rtol=1e-5))
    assert(np.allclose(spr.reflectance_ort, ref_R_ort, rtol=1e-5))


def test_mstmv4_halfspace():
    spr = SPR_v4(_in_periodic['wls'], temp_dir='./temp/')
    spr.environment_material = _in_periodic['n_env']
    D = _in_periodic['D']
    spheres = ExplicitSpheres(2, [-1.05*D, 0.0, 0.0, D/2.,
                                   1.05*D, 0.0, 0.0, D/2.],
                              mat_filename=2*[Material(_in_periodic['n'])])
    spr.set_spheres(spheres)
    spr.set_incident_field(fixed=True)
    spr.set_boundary(False)
    spr.set_layers([Material(1.0)])
    assert(not spr._check_spheres_in_layer())
    # move spheres deeper
    spheres = ExplicitSpheres(2, [-1.05*D, 0.0, -D/1.95, D/2.,
                                   1.05*D, 0.0, -D/1.95, D/2.],
                              mat_filename=2*[Material(_in_periodic['n'])])
    spr.set_spheres(spheres)
    assert(spr._check_spheres_in_layer())

    spr.simulate()

    import matplotlib.pyplot as plt
    plt.plot(spr.wavelengths, spr.extinction_up, label='ext_up')
    plt.plot(spr.wavelengths, spr.extinction_dn, label='ext_dn')
    plt.plot(spr.wavelengths, spr.extinction_up_par, label='ext_up_par')
    plt.plot(spr.wavelengths, spr.extinction_dn_par, label='ext_dn_par')
    plt.legend()
    plt.show()

    plt.plot(spr.wavelengths, spr.scattering_up, label='sca_up')
    plt.plot(spr.wavelengths, spr.scattering_dn, label='sca_dn')
    plt.plot(spr.wavelengths, spr.scattering_up_par, label='sca_up_par')
    plt.plot(spr.wavelengths, spr.scattering_dn_par, label='sca_dn_par')
    plt.legend()
    plt.show()





if __name__ == '__main__':
    # ~ test_mstmv4_periodicity()
    test_mstmv4_halfspace()

