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
 Start new fitting.
 Inital positions correspond to cubic lattice of particles
"""

import sys
# set the path to mstm_spectrum scripts. Binary mstm files should be in current folder.
sys.path.append('/home/leon/ltg_projects/fit-T-matrix/mstm-spectrum')

from mstm_spectrum import SPR, ExplicitSpheres, FilmBackground, LorentzBackground
import fit_spheres_optic
from fit_spheres_optic import *

data = read_ascii('optic_sample22.dat', True, 0) # read and sort by 0th column

if max(data[0,:]) < 10:
    print('Data X column is in microns, will rescale to nm.')
    data[0,:] = data[0,:] * 1000; # to nanometers

#wavelengths, exp = rebin(400, 800, 51, data[0,:], data[1,:])
wavelengths, exp = rebin(300, 800, 51, data[0,:], data[1,:])

fit_spheres_optic.MATRIX_MATERIAL = 'Glass'  # 1.66
#~ prepare_fit(wavelengths, exp, LorentzBackground(wavelengths)  )
prepare_fit(wavelengths, exp )

### SET INITIAL VALUES ###
values = [0.012, 0.01] # scale, bkg
#~ values = [0.012, 0.01, 0.01, 20] # scale, bkg = a1+a2*film(depth=a3)
#~ values = [0.012, 0.01, 0.01, 100] # scale, bkg = a1+Lorentz(a2,a3)
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
N = (len(values)-2)/4
print 'Number of spheres: ', N
print 'Number of degrees of freedom: ', len(values)
print 'Matrix material : %s' % fit_spheres_optic.MATRIX_MATERIAL
raw_input('Press enter')

### OPTIMIZE (FIT) VALUES ###
result = so.fmin( func=target_func, x0=values, callback=cbplot, xtol=0.0001, ftol=0.001, maxiter=1000, full_output=True, disp=True, args=(wavelengths, exp) )

### DEAL WITH RESULTS ###
print(result)
values = result[0]
print values

#~ y_fit = get_spectrum( wavelengths, values )
#~ plt.plot( wavelengths, exp, wavelengths, y_fit )
#~ #plt.axis([0, 1, 1.1*np.amin(s), 2*np.amax(s)])
#~ plt.xlabel('Wavelength, nm')
#~ plt.ylabel('Exctinction, a.u.')
#~ plt.show()

#~ print('Gold film depth: %.2f'%values[3])

raw_input('press enter')
