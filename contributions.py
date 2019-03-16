# -*- coding: utf-8 -*-
#
# ----------------------------------------------------- #
#                                                       #
#  This code is a part of T-matrix fitting project      #
#  Contributors:                                        #
#   L. Avakyan <laavakyan@sfedu.ru>                     #
#   K. Yablunovskiy <kirill-yablunovskii@mail.ru>       #
#                                                       #
# ----------------------------------------------------- #
"""
Contributions to UV/vis extinction spectra other
then obtained from MSTM.
"""
#TODO: add backgtound contributions and adjust corresponding GUI panels. Remove background from mstm_spectrum.py
#TODO: add Mie contribution
#TDOD: add Kirill's neural network contribution
from __future__ import print_function
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

# use input in both python2 and python3
try:
   input = raw_input
except NameError:
   pass
# use xrange in both python2 and python3
try:
    xrange
except NameError:
    xrange = range

try:
    from film_exctinction import gold_film_ex  # for gold film background
except:
    pass

try:
    from mie_theory import calculate_mie_spectra
except:
    pass


class Contribution(object):
    """
    Abstract class to account for contributions other then
    calculated by MSTM. Here should come all lightweight calculated
    contribtions: constant background, lorentz and guass peaks, Mie, etc.
    Each should be implemented as a child class.
    """
    number_of_params = 0  # Should be another value in child class

    def __init__(self, wavelengths=[], name='ExtraContrib'):
        """
        Parameter object

        name : str
        wavelengths : list or np.array
            wavelengths
        values : list
            adjustable parameters, like constant bakcground, or peak parameters
        """
        self.name = name
        self.set_wavelengths(wavelengths)

    def set_wavelengths(self, wavelengths):
        self.wavelengths = np.array(wavelengths)

    def calculate(self, values):
        """
        return np.array of calculated contribution at wavelength.

        This method should be overriden in child classes.
        """
        self._check(values)
        return np.zeros(len(self.wavelengths))

    def _check(self, values):
        if len(values) < self.number_of_params:
            raise Exception('Too few values! '+str(values))

    def plot(self, values, fig=None, axs=None):
        """
        plot contribution
        values : list of parameters
        fig, axs : matplotlib objects
        """
        flag = fig is None
        if flag:
            fig = plt.figure()
            axs = fig.add_subplot(111)
        x = self.wavelengths
        y = self.calculate(values)
        axs.plot(x, y, 'g--', label=self.name)
        axs.set_ylabel('Intensity')
        axs.set_xlabel('Wavelength, nm')
        axs.legend()
        if flag:
            plt.show()


class ConstantBackground(Contribution):
    """
    Simple background contribution ruled by single parameter.
    """
    number_of_params = 1

    def calculate(self, values):
        self._check(values)
        return values[0] * np.ones(len(self.wavelengths))


class LinearBackground(Contribution):
    """
    Two-parameter background $ a * {\lambda} + b $.
    """
    number_of_params = 2

    def calculate(self, values):
        self._check(values)
        return values[0] + values[1] * self.wavelengths


class LorentzPeak(Contribution):
    """
    Lorentz function
    """
    number_of_params = 3

    def calculate(self, values):
        self._check(values)
        return values[0] / ((self.wavelengths - values[1])**2 + (values[2])**2)


class GaussPeak(Contribution):
    """
    Gauss function
    """
    number_of_params = 3

    def calculate(self, values):
        self._check(values)
        return values[0] * np.exp(-(self.wavelengths - values[1])**2 / (2 * values[2]**2))


class LorentzBackground(Contribution):
    """
    Lorentz peak in background. Peak center is fixed.
    """
    number_of_params = 3

    def calculate(self, values, center=250):
        self._check(values)
        return values[0] + values[1] / ((self.wavelengths-center)**2 + (values[2])**2)


class FilmBackground(Contribution):
    """
        Background interpolated from experimental spectra of gold foil
    """
    number_of_params = 3

    def calculate(self, values):
        self._check(values)
        return values[0] + values[1] * gold_film_ex(values[2], self.wavelengths)


class MieSingleSphere(Contribution):
    """
        Mie contribution from single sphere
    """
    number_of_params = 2
    material = None  # instance of mstm_spectrum.Material
    matrix = 1.0     # medium refractive index

    def calculate(self, values):
        self._check(values)
        if self.material is None:
            raise Exception('Mie calculation requires material data. Stop.')
        _, _, mie_extinction, _ = calculate_mie_spectra(
            self.wavelengths, np.abs(values[1]), self.material, self.matrix
        )
        return values[0] * mie_extinction

    def set_material(self, material, matrix=1.0):
        self.matrix = matrix
        try:
            material.get_n(550)
            material.get_k(550)
        except:
            raise Exception('Bad material object')
        self.material = material


class MieLognormSpheres(MieSingleSphere):
    """
        Mie contribution from an ensemble of spheres
        with sizes distributed by Lognormal law
    """
    number_of_params = 3
    diameters = np.logspace(0, 3, 301)
    MAX_DIAMETER_TO_PLOT = 100

    def lognorm(self, x, mu, sigma):
        return (1.0/(x*sigma*np.sqrt(2*np.pi)))*np.exp(-((np.log(x)-mu)**2)/(2*sigma**2))

    def calculate(self, values):
        self._check(values)
        dD = np.ediff1d(self.diameters, to_begin=1e-3)
        distrib = self.lognorm(self.diameters, values[1], values[2])
        result = np.zeros_like(self.wavelengths)
        for diameter, count in zip(self.diameters, distrib*dD):
            self.number_of_params = 2  # else will get error on a check
            mie_ext = super(MieLognormSpheres, self).calculate(values=[1.0, diameter/2.0])
            self.number_of_params = 3  # ugly, but everything has a price
            result += count * mie_ext
        return values[0] * result

    def plot_distrib(self, values, fig=None, axs=None):
        """
        plot size distribution
        values : list of parameters
        fig, axs : matplotlib objects
        """
        flag = fig is None
        if flag:
            fig = plt.figure()
            axs = fig.add_subplot(111)
        x = self.diameters[self.diameters < self.MAX_DIAMETER_TO_PLOT]
        y = self.lognorm(x, values[1], values[2])
        axs.plot(x, y, 'b', label=self.name)
        axs.set_ylabel('Count')
        axs.set_xlabel('Diameter, nm')
        axs.legend()
        if flag:
            plt.show()


if __name__=='__main__':
    # tests come here
    cb = ConstantBackground(name='const', wavelengths=[300,400,500,600,700,800])
    print(cb.calculate([3]))
    cb.plot([3])
    print('See you!')

