#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#-----------------------------------------------------#
#                                                     #
# This code is a part of T-matrix fitting project     #
# Contributors:                                       #
#  L. Avakyan <laavakyan@sfedu.ru>,                   #
#  A. Skidanenko <ann.skidanenko@yandex.ru>           #
#                                                     #
#-----------------------------------------------------#
'''
  Example file to perform a fitting from scratch.
  Inital positions correspond to cubic lattice of particles
'''

import sys
# set the path to mstm_spectrum scripts. Binary mstm files should be in current folder.
sys.path.append('/home/leon/ltg_projects/fit-T-matrix/mstm-spectrum')

from mstm_spectrum import ExplicitSpheres
from alloy_AuAg import AlloyAuAg
from fit_spheres_optic import Fitter, FixConstraint



fitter = Fitter('example/optic_sample19.dat')
fitter.set_background('lorentz')
fitter.set_matrix('glass')

spheres = ExplicitSpheres(3, [-10, 10, 0], [0, 0, 14], [0, 0, 0],
                          [12, 12, 12], [AlloyAuAg(x_Au=0.0)]*3)
fitter.set_spheres(spheres)

fitter.add_constraint(FixConstraint('x0'))
fitter.add_constraint(FixConstraint('y0'))
fitter.add_constraint(FixConstraint('z0'))

fitter.report_freedom()
raw_input('Press enter')

fitter.run()

fitter.report_result()
raw_input('press enter')
