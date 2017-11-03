#!/usr/bin/python
# -*- coding: utf-8 -*-
#-----------------------------------------------------#
#                                                     #
# This code is a part of T-matrix fitting project     #
# Contributors:                                       #
#  L. Avakyan <laavakyan@sfedu.ru>                    #
#  A. Skidanenko <ann.skidanenko@ya.ru>               #
#                                                     #
#-----------------------------------------------------#
"""
  Fitting of particle aggreagate T-matrix spectrum
  to experimental SPR spectrum.
"""

import os
from mstm_spectrum import SPR, ExplicitSpheres, Background, FilmBackground
import numpy as np
from scipy import interpolate
import scipy.optimize as so

import matplotlib.pyplot as plt

class Parameter(object):
    """
    Class for parameter object used for storage of
    parameter's name, value and variation limits.
    """
    def __init__(self, name, value=1, min=None, max=None, internal_loop=False):
        """
        Parameter object

        name : str
            name of parameter used for constraints etc
        value : float
            initial value of parameter
        min, max : float
            bounds for parameter variation (optional)
        internal_loop : bool
            if True the parameter will be allowed to vary in internal (fast)
            loop, which does not require recalculation of spectrum
        """
        self.name = name
        self.value = self.ini_value = value
        self.min = min
        self.max = max
        # ...

    def __str__(self):
        return self.name

class Fitter(object):
    """
    Class to perform fit of experimental Exctinction spectrum
    """
    def __init__(self, exp_filename, wl_min=300, wl_max=800, wl_npoints = 51,
                 bkg_method='constant', plot_progress=True):
        """
        Creates the Fitter object

        Parameters:
        exp_filename : str
            name of file with experimental data
        wl_min, wl_max : float (nm)
            wavelength bounds for fitting.
        wl_npoints : int
            number of wavelengths where spectra will be calcualted and compared.
        plot_progress : bool
            show fitting progress using matplotlib.
            Should be turned off when run on parallel cluster.
        """
        self.exp_filename = exp_filename
        data = np.loadtxt(self.exp_filename)    # load data
        data = data[np.argsort(data[:,0]),:]    # sort by 0th column
        if np.max(data[:,0]) < 10:      # if values are really low
            print('WARINING: Data X column is probably in mum, automatilcally rescaling to nm.')
            data[:,0] = data[:,0] * 1000

        self.wl_min = max(wl_min, data[ 0, 0])
        self.wl_max = min(wl_max, data[-1, 0])
        self.wl_npoints = wl_npoints
        print('Wavelength limits are setted to: %f < wl < %f'% (self.wl_min, self.wl_max))
        self.wls, self.exp = self._rebin(self.wl_min, self.wl_max, self.wl_npoints,
                                         data[:,0], data[:,1])
        self.params = {}        # dictionaty of parameters objects
        self.calc = []          # calculated spectrum
        self.chidq = -1         # chi square (squared residual)

        self.plot_progress = plot_progress
        if self.plot_progress:
            plt.ion()
            self.fig = plt.figure()
            ax = self.fig.add_subplot(111)
            ax.plot(self.wls, self.exp, 'ro')
            self.line1, = ax.plot(self.wls, np.zeros_like(self.wls), 'b-')
            self.fig.canvas.draw()
        # set scale as default
        set.set_scale()
        # set background method as default
        self.background = None
        self.set_background()
        # set matrix material as default
        self.set_matrix()

    def _rebin(self, xmin, xmax, N, x, y):
        """
        hidden method used to rebin data to uniform scale
        """
        f = interpolate.interp1d(x, y)
        xnew = np.linspace(xmin, xmax, N)
        ynew = f(xnew)
        return xnew, ynew

    def set_matrix(self, material='AIR'):
        """
        set refraction index of matrix material

        material : {'AIR'|'WATER'|'GLASS'} or float
            the name of material or
            refraction index value.
        """
        self.MATRIX_MATERIAL = material

    def set_scale(self, value=0.1):
        if 'scale' is not in self.params.keys:
            self.params['scale'] = Parameter('scale', value=value, internal_loop=True)
        else:
            self.params['scale'].value = value
            self.params['scale'].ini_value = value

    def set_background(self, bkg_method='constant', initial_values = None):
        """
            Set method of background treatment.

            bkg_method : {'constant'|'linear'|'lorentz'|'gold_film'}
                Name of the method. If not valid than simple
                constant background will be used.
            initial_values : float array
        """
        # remove old bkg parameters
        if self.background is not None:
            n = self.background.number_of_params()
            for i in range(n):
                self.params.pop('bkg%i' % i)
        # create object
        bkg_method = bkg_method.lower()
        if bkg_method == 'linear':
            self.background = LinearBackground(self.wls)
        elif bkg_method == 'lorentz':
            self.background = LorentzBackground(self.wls)
        elif bkg_method == 'gold_film':
            self.background = FilmBackground(self.wls)
        else: # 'constant'
            self.background = Background(self.wls)
        print('Background method: %s' % self.background)
        # create new parameter objects
        n = self.background.number_of_params()
        for i in range(n):
            self.params['bkg%i' % i] = Parameter('bkg%i' % i, value=1, internal_loop=True)
            if initial_values is not None:
                assert len(initial_values) == n
                self.params['bkg%i' % i].value = initial_values[i]
                self.params['bkg%i' % i].ini_value = initial_values[i]

    def set_spheres(self, spheres):
        #TODO: implement it
        pass

    def get_spectrum(self, values):
        """
        Calculate the spectrum of agglomerates using mstm_spectrum module

        values : float array
            the values array passed from scipy optimizer
        """
        #print('Current scale: %f bkg: %f' % (values[0], values[1]))
        spr = SPR( self.wls )
        spr.environment_material = self.MATRIX_MATERIAL

        result = np.zeros(len(wavelengths))
        try:
            N = ( len(values)-1-background.number_of_params() )/4
            if N > 0:
                #print 'N=', N
                spr.set_spheres(  ExplicitSpheres( N, np.array(values[-N*4:]), [], [], [], mat_filename='etaGold.txt')  )
                _, extinction = spr.simulate('exct.dat')
                result = np.array(extinction)
            else:
                result = np.zeros(len(wavelengths))
        except Exception as e:
            print e
            result = np.zeros(len(wavelengths))
        calculated_extinction = values[0]*result + background.get_bkg(values[1:]) # scale * extinction + bkg
        return calculated_extinction

    def target_func(self, values, x_dat, y_dat):
        print( 'Scale: %f Bkg: %f'% (values[0], values[1]) )
        y_fit =  get_spectrum( wavelengths, values, MATRIX_MATERIAL )
        global chisq
        #chisq = np.sum( (y_fit**3 - y_dat**3)**2 / y_dat**3 )
        chisq = np.sum( (y_fit - y_dat)**2 * y_dat**3 ) * 1E3
        #print( chisq )
        return chisq

    def _cbplot(self, values ):
        """
        callback function
        """
        print( 'Scale: %0.3f Bkg: %0.3f ChiSq: %.8f'% (values[0], values[1], chisq) )
        if self.plot_progress:
            line1.set_ydata(calculated_extinction)
            fig.canvas.draw()


if __name__ == '__main__':
    print(Parameter('test_prm'))
    fitter = Fitter('example/optic_sample22.dat')

    #~ ### SET INITIAL VALUES ###
    #~ values = [0.02, 0.01] # scale, bkg
    #~ A = 200 # 400
    #~ a = 20
    #~ d = 100
    #~ x = -(A/2.0)
    #~ while x < (A/2.0):
        #~ y = -(A/2.0)
        #~ while y < (A/2.0):
            #~ z = -(A/2.0)
            #~ while z < (A/2.0):
                #~ if (x*x+y*y+z*z < A*A/4.0):
                    #~ values.append(x)
                    #~ values.append(y)
                    #~ values.append(z)
                    #~ values.append(a)
                    #~ #print x, y, z
                #~ z = z + (2*a+d)
            #~ y = y + (2*a+d)
        #~ x = x + (2*a+d)
    #~ N = (len(values)-2)/4
    #~ print 'Number of spheres: ', N
    #~ print 'Number of degrees of freedom: ', len(values)
    #~ MATRIX_MATERIAL = 'Glass'
    #~ print 'Matrix material : ', MATRIX_MATERIAL
    #~ raw_input('Press enter')
#~
    #~ ### OPTIMIZE (FIT) VALUES ###
    #~ result = so.fmin( func=target_func, x0=values, callback=cbplot, xtol=0.0001, ftol=0.001, maxiter=1000, full_output=True, disp=True, args=(wavelengths, exp) )
#~
    #~ ### DEAL WITH RESULTS ###
    #~ print(result)
    #~ values = result[0]
    #~ print values

    #~ y_fit = get_spectrum( wavelengths, values )
    #~ plt.plot( wavelengths, exp, wavelengths, y_fit )
    #~ #plt.axis([0, 1, 1.1*np.amin(s), 2*np.amax(s)])
    #~ plt.xlabel('Wavelength, nm')
    #~ plt.ylabel('Exctinction, a.u.')
    #~ plt.show()
    raw_input('press enter')
    print 'It is over.'
