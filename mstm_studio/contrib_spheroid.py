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
'''
Contributions to optical extinction spectra from axial-symmetric
particles. Currently, spheroids.
'''
import numpy as np
try:
    from scatterpy.tmatrix import calc_T
    from scatterpy.shapes import spheroid, chebyshev, gen_chebyshev
except ImportError:
    print('WARNING: Could not load `scatterpy` library!\n'
          'Spheroid functional will be disabled')

from mstm_studio.contributions import MieSingleSphere


class SpheroidSP(MieSingleSphere):
    '''
    Extinction from spheroid calculated in T-matrix approach
    using external library `ScatterPy`
    <https://github.com/TCvanLeth/ScatterPy>
    '''

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.number_of_params = 3
        self.NORDER = 5   # number of harmonics
        self.NGAUSS = 11  # integration points

    def _get_sfunc(self, values):
        return spheroid(np.array([np.abs(values[2])]))

    def calculate(self, values):
        '''
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
        '''
        self._check(values)
        if self.material is None:
            raise Exception('T-matrix calculation requires material data. Stop.')
        Cext = np.zeros(len(self.wavelengths))
        if calc_T is None:  # failed to import scatterpy
            return Cext
        self.material.D = np.abs(values[1])

        nk = self.material.get_nk(self.wavelengths)
        for iwl, wl in enumerate(self.wavelengths):
            size_param = 2 * np.abs(values[1] / 2.0) * self.matrix
            T = calc_T(size_param,
                       wl, nk[iwl] / self.matrix,  # rtol=0.001,
                       n_maxorder=self.NORDER, n_gauss=self.NGAUSS,
                       sfunc=lambda x: self._get_sfunc(values))
            Nmax = T.shape[-3]
            for n in range(1, Nmax+1):
                Cext[iwl] += np.real(T[0, 0, n-1, n-1, 0, 0] +
                                     T[0, 0, n-1, n-1, 1, 1])
                for m in range(1, n+1):
                    Cext[iwl] += 2 * np.real(T[0, m, n-1, n-1, 0, 0] +
                                             T[0, m, n-1, n-1, 1, 1])
        Cext = -self.wavelengths**2 / (2 * np.pi) * Cext
        Cext = Cext / (np.pi * size_param**2 / 4.0)
        Cext = Cext * self.matrix  # to compare with mstm's results
        return values[0] * Cext

    def plot_shape(self, values, fig=None, axs=None):
        '''
        Plot shape profile.
        Spatial shape is achieved by rotation over vertical axis.

        Parameters:

            values: list of control parameters
                    `scale`, `size` and `aspect`

            fig: matplotlib figure

            axs: matplotlib axes

        Return:

            filled/created fig and axs objects
        '''
        flag = fig is None
        if flag:
            fig = plt.figure()
            axs = fig.add_subplot(111)
        theta = np.linspace(0, 2*np.pi, 100)
        r, _ = self._get_sfunc(values)(np.cos(theta))
        r = np.squeeze(r)  # remove extra dimension
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


class ChebyshevSP(SpheroidSP):
    '''
    Extinction from Chebyshev shaped
    single particle calculated
    using external library `ScatterPy`
    <https://github.com/TCvanLeth/ScatterPy>

    Parameters:

        values: list of parameters `scale`, `size`,
                Chebyshev polynom order `n` and
                deformation parameter `eps`

        Scale is an arbitrary multiplier.

        Size parameter is the radius of undeformed
        sphere.

        Polynom order, positive int number.
        Sphere: n = 0,
        n = 1 - spheroid, etc

        Deformation parameter can be both positive and negative
        Should be small (<< 1) to comply with Raylaigh criterion (?).
        Try bigger values with tweaked NORDER and NGAUSS.

    '''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.number_of_params = 4  # scale, size, n, eps
        self.NORDER = 5   # number of harmonics
        self.NGAUSS = 15  # integration points

    def _get_sfunc(self, values):
        # values[0] -- multiplier
        # values[1] -- size
        n = int(np.round(values[2]))  # poly order
        # print(f'poly order: {n}')
        eps = values[3]  # deformation
        # TODO: auto increase NORDER for big eps ?
        return chebyshev(np.array([eps]), n)


if __name__ == '__main__':
    from mstm_studio.mstm_spectrum import Material
    import matplotlib.pyplot as plt
    import os

    n_env = 1.33
    mat_gold = Material(os.path.join('nk', 'etaGold.txt'))
    wls = np.linspace(300, 800, 45)
    npsize = 10  # diameter of nanoparticle

    if False:
        sph = SpheroidSP(wavelengths=wls)
        values = [1, npsize, 2]
    else:
        sph = ChebyshevSP(wavelengths=wls)
        values = [1, npsize, 2, 0.5]
        sph.NORDER = 30
        sph.NGAUSS = 30

    sph.set_material(mat_gold, n_env)
    sph.plot_shape(values)

    ext_sph = sph.calculate(values)
    sph.plot(values)



