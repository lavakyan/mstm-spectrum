# -*- coding: utf-8 -*-
#
# ----------------------------------------------------- #
#                                                       #
#  This code is a part of T-matrix fitting project      #
#  Contributors:                                        #
#   L. Avakyan <laavakyan@sfedu.ru>                     #
#   D. Kostyulin <kostyulin@sfedu.ru>                   #
#   E. Mamyan <emamian@sfedu.ru>                        #
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
        self.nm_max = 100  # for autotuning, if NORDER is none
        self.ng_max = 500  # for autotuning, if NGAUSS is none

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
                       nm_max=self.nm_max, ng_max=self.ng_max,
                       sfunc=lambda x: self._get_sfunc(values))
            Nmax = T.shape[-3]
            if T.ndim == 6:
                T = T[0]
            for n in range(1, Nmax+1):
                Cext[iwl] += np.real(T[0, n-1, n-1, 0, 0] +
                                     T[0, n-1, n-1, 1, 1])
                for m in range(1, n+1):
                    Cext[iwl] += 2 * np.real(T[m, n-1, n-1, 0, 0] +
                                             T[m, n-1, n-1, 1, 1])
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

    Chebyshev shaped particles are rotationally (phi)
    symmetric particles with angular dependence:

    r(θ) = r₀ [1 + eps Tₙ(cos θ)]

    where Tₙ(cos θ) is a Chebyshev polynom,
    Tₙ(cos θ) = cos (n θ)

    Parameters:

        values: list of parameters `scale`, `size`,
                Chebyshev polynom order `n` and
                deformation parameter `eps`

        Scale is an arbitrary multiplier.

        Size parameter is the radius of undeformed
        sphere.

        Polynom order, positive int number.
        Sphere: n = 0, eps ignored
        n = 2 - spheroid, etc

        Deformation parameter eps can be both positive and negative
        and |eps| < 1.
        Should be small (<< 1) for better convergence.
        Check results with respect to NORDER and NGAUSS
        for strong deformations!

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
        if eps < -0.999:
            print(f'WARNING: to big deformation! {eps}')
            eps = -0.999
        elif eps > 0.999:
            eps = 0.999
            print(f'WARNING: to big deformation! {eps}')
        # TODO: auto increase NORDER for big eps ?
        return chebyshev(np.array([eps]), n)


class GenChebyshevSP(SpheroidSP):
    '''
    Extinction from "generalized" Chebyshev shaped
    single particle calculated
    using external library `ScatterPy`
    <https://github.com/TCvanLeth/ScatterPy>

    Chebyshev shaped particles are rotationally (phi)
    symmetric particles with angular dependence:

    r(θ) = r₀ sum_n [1 + eps Tₙ(cos θ)]

    where Tₙ(cos θ) is a Chebyshev polynom,
    Tₙ(cos θ) = cos (n θ)

    Parameters:

        values: list of parameters `scale`, `size`,
                deformation parameters

        Scale is an arbitrary multiplier.

        Size parameter is the radius of undeformed
        sphere.

        Polynom order is dedicted from the length of
        deformation parameters list, i.e. as
        `len(values) - 2 + 1`

        Hard-coded limitation of 10 order.

        Any deformation parameter can be both positive and negative
        and limited by |eps| < 1.
        They should be small (<< 1) for better convergence.
        Check results with respect to NORDER and NGAUSS
        for strong deformations!
    '''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.number_of_params = 1 + 1 + 10  # scale, size, deformations
        self.NORDER = 5   # number of harmonics
        self.NGAUSS = 15  # integration points

    def _get_sfunc(self, values):
        # values[0] -- multiplier
        # values[1] -- size
        epss = np.array(values[2:])  # list of deformations
        n = len(epss)
        # print(f'poly order: {n}')
        # print(epss)
        epss[epss < -0.999] = -0.999
        epss[epss > 0.999] = 0.999
        # TODO: auto increase NORDER for big eps ?
        return gen_chebyshev(epss, ng=60)

    def calculate(self, values):
        _values = np.zeros(self.number_of_params)
        _values[:len(values)] = values
        # print(_values)
        return super().calculate(_values)


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
    elif False:
        sph = ChebyshevSP(wavelengths=wls)
        values = [1, npsize, 5, 0.15]
        sph.NORDER = 30
        sph.NGAUSS = 30
    else:
        sph = GenChebyshevSP(wavelengths=wls)
        values = [1, npsize, 0.0, 0.0, 0.0, 0.2, -0.1]
        sph.NORDER = 50
        sph.NGAUSS = 50

    sph.set_material(mat_gold, n_env)
    sph.plot_shape(values)

    ext_sph = sph.calculate(values)
    sph.plot(values)



