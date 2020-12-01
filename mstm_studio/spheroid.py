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
            print('SpheroidSP: current wavelength %.0f nm' % wl)
            T = calc_T(values[1], wl, nk[iwl] / self.matrix,  # rtol=0.001,
                       n_maxorder=7, n_gauss=7*2+1,
                       sfunc=lambda x: spheroid(np.array([values[2]])))
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
    sph = SpheroidSP(wavelengths=np.linspace(300, 800, 51))
    sph.set_material(AlloyAuAg(x_Au=1), 1.6)
    sph.plot([1, 10, 1.0])  # scale diameter aspect

