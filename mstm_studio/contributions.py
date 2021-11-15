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
from __future__ import print_function
from __future__ import division
import numpy as np
try:
    import matplotlib.pyplot as plt
except:
    pass

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
    from mstm_studio.mie_theory import calculate_mie_spectra
except:
    pass


class Contribution(object):
    """
    Abstract class to include contributions other then
    calculated by MSTM. All lightweight calculated
    contribtions (constant background, lorentz and guass peaks, Mie, etc.)
    should enhirit from it.
    """
    number_of_params = 0  # Should be another value in child class

    def __init__(self, wavelengths=[], name='ExtraContrib'):
        """
        Parameters:

            wavelengths: list or numpy array
                wavelengths in nm

            name: string
                optional label

        """
        self.name = name
        self.set_wavelengths(wavelengths)

    def set_wavelengths(self, wavelengths):
        """
        Modify wavelengths
        """
        self.wavelengths = np.array(wavelengths)

    def calculate(self, values):
        """
        This method should be overriden in child classes.

        Parameters:

            values: list of control parameters

        Return:

            numpy array of contribution values at specified wavelengths
        """
        self._check(values)
        return np.zeros(len(self.wavelengths))

    def _check(self, values):
        if len(values) < self.number_of_params:
            raise Exception('Too few values! '+str(values))

    def plot(self, values, fig=None, axs=None):
        """
        plot contribution

        Parameters:

            values: list of parameters

            fig: matplotlib figure

            axs: matplotlib axes

        Return:

            filled/created fig and axs objects
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
        return fig, axs


class ConstantBackground(Contribution):
    """
    Constant background contribution :math:`f(\lambda) = bkg`.
    """
    number_of_params = 1

    def calculate(self, values):
        """
        Parameters:

            values: [bkg]

        Return:

            numpy array
        """
        self._check(values)
        return values[0] * np.ones(len(self.wavelengths))


class LinearBackground(Contribution):
    """
    Two-parameter background :math:`f(\lambda) = a \cdot \lambda + b`.
    """
    number_of_params = 2

    def calculate(self, values):
        """
        Parameters:

            values: list of control parameters `scale`, `mu` and `Gamma`

        Return:

            numpy array
        """
        self._check(values)
        return values[0] + values[1] * self.wavelengths


class LorentzPeak(Contribution):
    """
    Lorentz function

    .. math::

        L(\lambda) = \\frac {scale} {(\lambda-\mu)^2 + \Gamma^2}

    """
    number_of_params = 3

    def calculate(self, values):
        """
        Parameters:

            values: list of control parameters `scale`, `mu` and `Gamma`

        Return:

            numpy array
        """
        self._check(values)
        return values[0] / ((self.wavelengths - values[1])**2 + (values[2])**2)


class GaussPeak(Contribution):
    """
    Gauss function

    .. math::

        G(\lambda) = scale \cdot \exp\left( - \\frac{(\lambda-\mu)^2}{2\sigma^2} \\right)

    """
    number_of_params = 3

    def calculate(self, values):
        """
        Parameters:

            values: list of control parameters `scale`, `mu` and `sigma`

        Return:

            numpy array
        """
        self._check(values)
        return values[0] * np.exp(-(self.wavelengths - values[1])**2 / (2 * values[2]**2))


class LorentzBackground(Contribution):
    """
    Lorentz peak in background. Peak center is fixed.

    .. math::

        L(\lambda) = \\frac {scale} {(\lambda-center)^2 + \Gamma^2}

    """
    number_of_params = 2
    center = 250

    def calculate(self, values):
        self._check(values)
        return values[0] / ((self.wavelengths-self.center)**2 + (values[1])**2)


class FilmBackground(Contribution):
    """
    Background interpolated from experimental spectra of gold foil
    """
    number_of_params = 3

    def calculate(self, values):
        """ TODO """
        self._check(values)
        return values[0] + values[1] * gold_film_ex(values[2], self.wavelengths)


class MieSingleSphere(Contribution):
    """
    Mie contribution from single sphere.

    Details are widely discusses, see, for example [Kreibig_book1995]_
    """
    number_of_params = 2
    material = None  # instance of mstm_spectrum.Material
    matrix = 1.0     # medium refractive index

    def calculate(self, values):
        """
        Parameters:

            values: list of control parameters `scale`, `diameter`

        Return:

            extinction efficiency array of Mie sphere
        """
        self._check(values)
        if self.material is None:
            raise Exception('Mie calculation requires material data. Stop.')
        D = np.abs(values[1])
        self.material.D = D
        _, _, mie_extinction, _ = calculate_mie_spectra(
            self.wavelengths, D, self.material, self.matrix)
        return values[0] * mie_extinction

    def set_material(self, material, matrix=1.0):
        """
        Define the material of sphere and environment

        Parameters:

            material: Material object
                material of the sphere

            matrix: float, string or Material object
                material of the environment

        Return:

            True if properties were changed, False - otherwise.

        """
        changed = False
        try:
            matr = float(matrix)
        except:
            matr = matrix.get_n(550)  # assume it is Material instance
        if matr != self.matrix:
            self.matrix = matr
            changed = True
        try:
            material.get_n(550)
            material.get_k(550)
        except:
            raise Exception('Bad material object')
        if material is not self.material:
            self.material = material
            changed = True
        return changed


class MieLognormSpheres(MieSingleSphere):
    """
    Mie contribution from an ensemble of spheres
    with sizes distributed by Log-Normal law
    """
    number_of_params = 3
    diameters = np.logspace(0, 3, 301)
    MAX_DIAMETER_TO_PLOT = 100

    def lognorm(self, x, mu, sigma):
        """
        The shape of Log-Normal distribution:

        .. math::

            LN(D) = \\frac {1}{D \sigma \sqrt{2\pi}} \exp\left( - \\frac{(\log(D)-\mu)^2}{2\sigma^2} \\right)
        """
        return (1.0/(x*sigma*np.sqrt(2*np.pi)))*np.exp(-((np.log(x)-mu)**2)/(2*sigma**2))

    def calculate(self, values):
        """
        Parameters:

            values: list of control parameters `scale`, `mu` and `sigma`

        Return:

            Mie extinction efficiency of log-normally distributed spheres
        """
        self._check(values)
        dD = np.ediff1d(self.diameters, to_begin=1e-3)
        distrib = self.lognorm(self.diameters, np.abs(values[1]), np.abs(values[2]))
        result = np.zeros_like(self.wavelengths)
        for diameter, count in zip(self.diameters, distrib*dD):
            self.number_of_params = 2  # else will get error on a check
            mie_ext = super(MieLognormSpheres, self).calculate(values=[1.0, diameter/2.0])
            self.number_of_params = 3  # ugly, but everything has a price
            result += count * mie_ext
        return values[0] * result

    def plot_distrib(self, values, fig=None, axs=None):
        """
        Plot size distribution

        Parameters:

            values: list of control parameters

            fig: matplotlib figure

            axs: matplotlib axes

        Return:

            filled/created fig and axs objects
        """
        flag = fig is None
        if flag:
            fig = plt.figure()
            axs = fig.add_subplot(111)
        x = self.diameters[self.diameters < self.MAX_DIAMETER_TO_PLOT]
        y = self.lognorm(x, np.abs(values[1]), np.abs(values[2]))
        axs.plot(x, y, 'b', label=self.name)
        axs.set_ylabel('Count')
        axs.set_xlabel('Diameter, nm')
        axs.legend()
        if flag:
            plt.show()
        return fig, axs

    def set_material(self, material, matrix=1.0):
        print('matrix: %s' % matrix)
        if super().set_material(material, matrix):
            self._M = None  # clear cache if changed materials


class MieLognormSpheresCached(MieLognormSpheres):
    """
    Mie contribution from an ensemble of spheres
    with sizes distributed by Lognormal law.

    Cached version - use it to speed-up fitting.
    """
    number_of_params = 3
    diameters = np.logspace(0, 3, 301)
    MAX_DIAMETER_TO_PLOT = 100
    _M = None  # cache matrix, None after initialization

    def lognorm(self, x, mu, sigma):
        return (1.0/(x*sigma*np.sqrt(2*np.pi)))*np.exp(-((np.log(x)-mu)**2)/(2*sigma**2))

    def calculate(self, values):
        """
        Parameters:

            values: list of control parameters `scale`, `mu` and `sigma`

        Return:

            Mie extinction efficiency of log-normally distributed spheres
        """
        self._check(values)
        if self._M is None:  # initialize cache matrix
            print('Building cache...')
            self._M = np.zeros(shape=(len(self.wavelengths), len(self.diameters)))
            self.number_of_params = 2  # else will get error on a check
            for i, diameter in enumerate(self.diameters):
                self._M[:, i] = super(MieLognormSpheres, self).calculate(values=[1.0, diameter/2.0])
            self.number_of_params = 3  # ugly, but everything has a price
            print('Building cache... done')

        dD = np.ediff1d(self.diameters, to_begin=1e-3)
        distrib = self.lognorm(self.diameters, np.abs(values[1]), np.abs(values[2]))
        result =  np.dot(self._M, distrib * self.diameters**2 * dD) / np.sum(distrib * self.diameters**2 * dD)
        return values[0] * result


if __name__=='__main__':
    # tests come here
    # ~ cb = ConstantBackground(name='const', wavelengths=[300,400,500,600,700,800])
    # ~ print(cb.calculate([3]))
    # ~ cb.plot([3])
    # ~ mie = MieSingleSphere(name='mie', wavelengths=np.linspace(300,800,50))
    # ~ mie = MieLognormSpheres(name='mie', wavelengths=np.linspace(300,800,50))
    from mstm_studio.contributions import MieLognormSpheresCached
    from mstm_studio.alloy_AuAg import AlloyAuAg
    mie = MieLognormSpheresCached(name='mie', wavelengths=np.linspace(300,800,50))
    mie.set_material(AlloyAuAg(x_Au=1), 1.66)
    mie.plot([1,1.5,0.5])  # scale mu sigma
    print('See you!')

