# -*- coding: utf-8 -*-
#
# ----------------------------------------------------- #
#                                                       #
#  This code is a part of T-matrix fitting project      #
#  Contributors:                                        #
#   L. Avakyan <laavakyan@sfedu.ru>                     #
#   D. Kostyulin <kostyulin@sfedu.ru>                   #
#                                                       #
# ----------------------------------------------------- #
"""
Contributions to optical extinction spectra from axial-symmetric
particles. Currently, spheroids.
"""
from __future__ import print_function
from __future__ import division
import numpy as np
try:
    import matplotlib.pyplot as plt
except ImportError:
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
    from scatterpy.tmatrix import calc_T  # , calc_T_inner
    from scatterpy.shapes import spheroid
except ImportError:
    print('WARNING: Could not load `scatterpy` library!')
    print('Spheroid functional will be disabled')
    pass

from mstm_studio.contributions import MieSingleSphere


class SpheroidSP(MieSingleSphere):
    """
    Extinction from spheroid calculated in T-matrix approach
    using external library `ScatterPy`
    <https://github.com/TCvanLeth/ScatterPy>
    """
    number_of_params = 3
    NORDER = 5   # number of harmonics
    NGAUSS = 11  # integration points

    def calculate(self, values):
        """
        Parameters:

            values: list of parameters `scale`, `size` and `aspect`
                    Scale is an arbitrary multiplier.
                    Size parameter is the radius of equivelent-volume
                    sphere.
                    The aspect ratio is
                    "the ratio of horizontal to rotational axes"
                    according to scatterpy/shapes.py

        Return:

            extinction efficiency array for spheroid particle
        """
        self._check(values)
        if self.material is None:
            raise Exception('T-matrix calculation requires material data. Stop.')
        Cext = np.zeros(len(self.wavelengths))
        if calc_T is None:  # failed to import scatterpy
            return Cext
        nk = self.material.get_n(self.wavelengths) + \
             1j * self.material.get_k(self.wavelengths)
        for iwl, wl in enumerate(self.wavelengths):
            # print('SpheroidSP: current wavelength %.0f nm' % wl)
            size_param = 2 * np.abs(values[1]) * self.matrix
            T = calc_T(size_param, wl, nk[iwl] / self.matrix,  # rtol=0.001,
                       n_maxorder=self.NORDER, n_gauss=self.NGAUSS,
                       sfunc=lambda x: spheroid(np.array([np.abs(values[2])])))
            Nmax = T.shape[-3]
            for n in range(1, Nmax+1):
                Cext[iwl] += np.real(T[0, 0, n-1, n-1, 0, 0] +
                                     T[0, 0, n-1, n-1, 1, 1])
                for m in range(1, n+1):
                    Cext[iwl] += 2 * np.real(T[0, m, n-1, n-1, 0, 0] +
                                             T[0, m, n-1, n-1, 1, 1])
        Cext = -self.wavelengths**2 / (2 * np.pi) * Cext
        Cext = Cext / (np.pi * size_param**2 / 4.0)
        return values[0] * Cext

    def plot_shape(self, values, fig=None, axs=None):
        """
        Plot shape profile.
        Spatial shape is achieved by rotation over vertical axis.

        Parameters:

            values: list of control parameters
                    `scale`, `size` and `aspect`

            fig: matplotlib figure

            axs: matplotlib axes

        Return:

            filled/created fig and axs objects
        """
        flag = fig is None
        if flag:
            fig = plt.figure()
            axs = fig.add_subplot(111)
        theta = np.linspace(0, 2*np.pi, 100)
        # 4/3 pi size^3 = 4/3 pi a^2 c
        # aspect = a / c
        a = values[1] * values[2]**(1/3.)
        c = a / values[2]
        # oblate / prolate speroid surface function from [Tsang1984]
        r = 1 / np.sqrt((np.sin(theta)/a)**2 + (np.cos(theta)/c)**2)
        x = r * np.sin(theta)
        z = r * np.cos(theta)
        axs.plot(x, z, 'b')
        axs.plot(-x, z, 'b')
        axs.plot([0, 0], [np.min(z), np.max(z)], 'b--')
        axs.set_aspect('equal', adjustable='box')
        axs.set_xlabel('X, nm')
        axs.set_ylabel('Z, nm')
        if flag:
            plt.show()
        return fig, axs

