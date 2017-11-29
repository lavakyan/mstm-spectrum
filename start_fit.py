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
"""
 Example file to perform a fitting from scratch.
 Inital positions correspond to cubic lattice of particles
"""

import sys
# set the path to mstm_spectrum scripts. Binary mstm files should be in current folder.
sys.path.append('/home/leon/ltg_projects/fit-T-matrix/mstm-spectrum')

from mstm_spectrum import ExplicitSpheres  #, FilmBackground, LorentzBackground
from fit_spheres_optic import Fitter


fitter = Fitter('example/optic_sample19.dat')
fitter.set_matrix('glass')


### SET INITIAL VALUES ###
values = []
A = 60   # 'box' size
a =  8   # sphere radius
d = 10   # 'gap' between spheres
x = -(A/2.0)
while x < (A/2.0):
    y = -(A/2.0)
    while y < (A/2.0):
        z = -(A/2.0)
        while z < (A/2.0):
            if (x*x+y*y+z*z < A*A/4.0):
                values.append(x)
                values.append(y)
                values.append(z)
                values.append(a)
                #print x, y, z
            z = z + (2*a+d)
        y = y + (2*a+d)
    x = x + (2*a+d)
N = len(values)/4

spheres = ExplicitSpheres(N, values, [], [], [], ['etaSilver.txt']*N)

fitter.set_spheres(spheres)

fitter.report_freedom()
raw_input('Press enter')

fitter.run()

fitter.report_result()
raw_input('press enter')
