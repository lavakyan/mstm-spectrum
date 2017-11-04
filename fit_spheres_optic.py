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
from mstm_spectrum import (SPR, ExplicitSpheres, Background,
        LinearBackground, FilmBackground, LorentzBackground)
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
        return '%s:%f' % (self.name, self.value)

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
        self.spheres = None
        self.calc = np.zeros_like(self.wls)  # calculated spectrum
        self.chidq = -1         # chi square (squared residual)

        self.plot_progress = plot_progress
        if self.plot_progress:
            plt.ion()
            self.fig = plt.figure()
            ax = self.fig.add_subplot(111)
            ax.plot(self.wls, self.exp, 'ro')
            self.line1, = ax.plot(self.wls, self.calc, 'b-')
            self.fig.canvas.draw()
        # set scale as default
        self.set_scale()
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
        if 'scale' not in self.params:
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
        """
        specify the spheres aggregate

        spheres : mstm_spectrum.Spheres
        """
        if self.spheres is not None:  # remove parameters of old spheres
            for i in xrange(self.spheres.N):
                self.params.pop('a%i' % i)
                self.params.pop('x%i' % i)
                self.params.pop('y%i' % i)
                self.params.pop('z%i' % i)

        self.spheres = spheres
        for i in xrange(self.spheres.N):
            self.params['a%i' % i] = Parameter('a%i' % i, self.spheres.a[i])
            self.params['x%i' % i] = Parameter('x%i' % i, self.spheres.x[i])
            self.params['y%i' % i] = Parameter('y%i' % i, self.spheres.y[i])
            self.params['z%i' % i] = Parameter('z%i' % i, self.spheres.z[i])

    def report_freedom(self):
        print('Number of spheres:\t%i' % len(self.spheres))
        Nbkg = self.background.number_of_params()
        print('Background parameters:\t%i' % Nbkg)
        #TODO: constraints!
        print('Degree of freedom:\t%i' % (1 + 4*N + Nbkg))  # scale + ..

    def update_spheres(self):
        """
        Set spheres radii and positions to values from params dict
        """
        assert self.spheres is not None
        for i in xrange(len(self.spheres)):
            self.spheres.a[i] = self.params['a%i' % i].value
            self.spheres.x[i] = self.params['x%i' % i].value
            self.spheres.y[i] = self.params['y%i' % i].value
            self.spheres.z[i] = self.params['z%i' % i].value

    def get_spectrum(self):
        """
        Calculate the spectrum of agglomerates using mstm_spectrum module.
        """
        #print('Current scale: %f bkg: %f' % (values[0], values[1]))
        #TODO: apply constraints on params
        spr = SPR(self.wls)
        spr.environment_material = self.MATRIX_MATERIAL
        result = np.zeros_like(self.wls)
        self.update_spheres()
        spr.set_spheres(self.spheres)
        try:
            _, extinction = spr.simulate('exct.dat')
            result = np.array(extinction)
        except Exception as e:
            print e
        bkg_values = []
        for key in self.params:
            if key.startswith('bkg'):
                bkg_values.append(self.params[key].value)
        self.calc = self.params['scale'].value * result + self.background.get_bkg(bkg_values)
        return self.calc

    def target_func(self, values):
        # unpack values to params
        i = 0
        for key in self.params:
            # check if varied here
            # check if not fixed, etc
            self.params[key].value = values[i]
            i += 1
        print('Scale: %f Bkg: %f' % (self.params['scale'].value, self.params['bkg0'].value) )
        y_dat = self.exp
        y_fit = self.get_spectrum()
        #~ self.chisq = np.sum( (y_fit**3 - y_dat**3)**2 / y_dat**3 )
        self.chisq = np.sum( (y_fit - y_dat)**2 )
        #~ self.chisq = np.sum( (y_fit - y_dat)**2 * y_dat**3 ) * 1E3
        #print(chisq)
        return self.chisq

    def _cbplot(self, values):
        """
        callback function
        """
        print('Scale: %0.3f Bkg: %0.3f ChiSq: %.8f'% (self.params['scale'].value,
              self.params['bkg0'].value, self.chisq) )
        if self.plot_progress:
            self.line1.set_ydata(self.calc)
            self.fig.canvas.draw()
            plt.pause(0.05)

    def run(self):
        # pack parameters to values
        values = []
        for key in self.params:
            #if belong to internal/external loop ..
            values.append(self.params[key].value)
        # run optimization
        result = so.fmin(func=self.target_func, x0=values, callback=self._cbplot, xtol=0.0001, ftol=0.001, maxiter=1000)
        ### DEAL WITH RESULTS ###
        #~ print(result)
        #~ values = result[0]
        #~ print values

if __name__ == '__main__':
    print(Parameter('test_prm'))
    fitter = Fitter('example/optic_sample22.dat')
    fitter.set_matrix('glass')
    #~ fitter.set_background('lorentz', [10, 20, 90])
    #~ ### SET INITIAL VALUES ###
    values = [] # coords and radii of spheres
    A = 200 # 400
    a = 20
    d = 100
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
    spheres = ExplicitSpheres(N, values)
    fitter.set_spheres(spheres)
    fitter.report_freedom()
    raw_input('Press enter')

    fitter.run()

    #~ y_fit = get_spectrum( wavelengths, values )
    #~ plt.plot( wavelengths, exp, wavelengths, y_fit )
    #~ #plt.axis([0, 1, 1.1*np.amin(s), 2*np.amax(s)])
    #~ plt.xlabel('Wavelength, nm')
    #~ plt.ylabel('Exctinction, a.u.')
    #~ plt.show()
    raw_input('press enter')
    print 'It is over.'
