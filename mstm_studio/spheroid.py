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
    from scatterpy.tmatrix import calc_T  #, calc_T_inner
    from scatterpy.shapes import spheroid
except:
    print('WARNING: Could not load `scatterpy` library!')
    print('Spheroid functional will be diabled')
    pass

from mstm_studio.contributions import MieSingleSphere


class SpheroidSP(MieSingleSphere):
    """
    Extinction from spheroid calculated in T-matrix approach
    using external library `ScatterPy`
    <https://github.com/TCvanLeth/ScatterPy>
    """
    number_of_params = 3
    NORDER = 7   # number of harmonics
    NGAUSS = 15  # integration points

    def calculate(self, values):
        """
        Parameters:

            values: list of parameters `scale`, `diameter` and `c`
                    The last is aspect ratio defined as
                    the ratio of horizontal to rotational axes.

        Return:

            extinction efficiency array for spheroid
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
            T = calc_T(size_param, wl, nk[iwl]/self.matrix,  # rtol=0.001,
                       n_maxorder=self.NORDER, n_gauss=self.NGAUSS,
                       sfunc=lambda x: spheroid(np.array([np.abs(values[2])])))
            Nmax = T.shape[-3]
            for n in range(1, Nmax+1):
                Cext[iwl] += np.real(T[0, 0, n-1, n-1, 0, 0] + \
                                     T[0, 0, n-1, n-1, 1, 1])
                for m in range(1, n+1):
                    Cext[iwl] += 2 * np.real(T[0, m, n-1, n-1, 0, 0] + \
                                             T[0, m, n-1, n-1, 1, 1])
        Cext = -self.wavelengths**2 / (2 * np.pi) * Cext
        return values[0] * Cext


if __name__=='__main__':
    # tests come here
    from mstm_studio.alloy_AuAg import AlloyAuAg
    wls = np.linspace(300, 800, 51)
    sph = SpheroidSP(wavelengths=wls)
    sph.set_material(AlloyAuAg(x_Au=1), 1.5)
    sph.NORDER = 5
    ext_sph = sph.calculate([1, 100, 1.0])
    #sph.plot([1, 10, 1.0])  # scale, diameter, aspect

    from mstm_studio.contributions import MieSingleSphere
    mie = MieSingleSphere(name='mie', wavelengths=wls)
    mie.set_material(AlloyAuAg(x_Au=1), 1.5)
    ext_mie = mie.calculate([1, 100])

    ext_sph = ext_sph / np.sum(ext_sph)
    ext_mie = ext_mie / np.sum(ext_mie)
    plt.plot(wls, ext_sph, label='T-mat')
    plt.plot(wls, ext_mie, label='Mie')
    plt.legend()
    plt.show()
