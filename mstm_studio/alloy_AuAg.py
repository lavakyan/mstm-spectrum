# -*- coding: utf-8 -*-
#
# ----------------------------------------------------- #
#                                                       #
#  This code is a part of T-matrix fitting project      #
#  Contributors:                                        #
#    A. Skidanenko <ann.skidanenko@ya.ru>               #
#    L. Avakyan <laavakyan@sfedu.ru>                    #
#    V. Pryadchenko <vvpryadchenko@sfedu.ru>            #
#                                                       #
# ----------------------------------------------------- #
"""
Dielectric function model for silver-gold alloy
according to Fortran code published in:
D. Rioux, S. Vallières, S. Besner, P. Muñoz, E. Mazur, and M. Meunier,
"An Analytic Model for the Dielectric Function of Au, Ag, and
 their Alloys" Adv. Opt. Mater. (2014) *2* 176-182
<http://dx.doi.org/10.1002/adom.201300457>
"""
from __future__ import print_function
import numpy as np
from mstm_studio.mstm_spectrum import Material
try:
    import matplotlib.pyplot as plt
except:
    pass

# use input in both python2 and python3
try:
    input = raw_input
except NameError:
    pass

c = 2.99792458e17
h = 4.135667516e-15
wp11 = 8.9234
wp12 = 8.5546
wp13 = 9.0218
gammap21 = 0.042389
gammap22 = 0.022427
gammap23 = 0.16713
einf31 = 2.2715
einf32 = 1.7381
einf33 = 2.2838
wg141 = 2.6652
wg142 = 4.0575
wg143 = 3.0209
w0151 = 2.3957
w0152 = 3.9260
w0153 = 2.7976
gamma161 = 0.17880
gamma162 = 0.017723
gamma163 = 0.18833
A171 = 73.251
A172 = 51.217
A173 = 22.996
w0281 = 3.5362
w0282 = 4.1655
w0283 = 3.3400
gamma291 = 0.35467
gamma292 = 0.18819
gamma293 = 0.68309
A2101 = 40.007
A2102 = 30.770
A2103 = 57.540


class AlloyAuAg(Material):
    """
    Material class for AuAg alloys.

    Use `get_n()` and `get_k()` to obtain values of refraction indexes (real
    and imaginary) at arbitraty wavelength (in nm) by model and code
    from Rioux et al doi:10.1002/adom.201300457
    """

    def __init__(self, x_Au):
        '''
        Parameters:

            x_Au: float
                fraction of gold
        '''
        self.x_Au = float(x_Au)
        self.__name__ = 'Mat_alloyAu%.2fAg%.2f' % (self.x_Au, 1-self.x_Au)

    def _drude(self, omega, GMF):
        wp = (GMF**2)*(2*wp11-4*wp13+2*wp12)+GMF*(-wp11+4*wp13-3*wp12)+wp12
        gammap = (GMF**2)*(2*gammap21-4*gammap23+2*gammap22)+GMF*(-gammap21+4*gammap23-3*gammap22)+gammap22
        einf = (GMF**2)*(2*einf31-4*einf33+2*einf32)+GMF*(-einf31+4*einf33-3*einf32)+einf32
        wg1 = (GMF**2)*(2*wg141-4*wg143+2*wg142)+GMF*(-wg141+4*wg143-3*wg142)+wg142
        w01 = (GMF**2)*(2*w0151-4*w0153+2*w0152)+GMF*(-w0151+4*w0153-3*w0152)+w0152
        gamma1 = (GMF**2)*(2*gamma161-4*gamma163+2*gamma162)+GMF*(-gamma161+4*gamma163-3*gamma162)+gamma162
        A1 = (GMF**2)*(2*A171-4*A173+2*A172)+GMF*(-A171+4*A173-3*A172)+A172
        w02 = (GMF**2)*(2*w0281-4*w0283+2*w0282)+GMF*(-w0281+4*w0283-3*w0282)+w0282
        gamma2 = (GMF**2)*(2*gamma291-4*gamma293+2*gamma292)+GMF*(-gamma291+4*gamma293-3*gamma292)+gamma292
        A2 = (GMF**2)*(2*A2101-4*A2103+2*A2102)+GMF*(-A2101+4*A2103-3*A2102)+A2102
        temp1 = omega**2 + 1j * omega*gammap
        drude = einf - ((wp**2)/temp1)
        compl2 = omega + 1j * gamma2
        compl1 = omega + 1j * gamma1
        cp1 = A1*((1/(compl1**2))*(-np.sqrt(compl1-wg1)*np.arctan(np.sqrt((wg1-w01)/(compl1-wg1)))) + (1/(compl1**2))*(-np.sqrt(compl1+wg1)*np.arctanh(np.sqrt((wg1-w01)/(compl1+wg1))))+(1/(compl1**2))*(2*np.sqrt(wg1)*np.arctanh(np.sqrt((wg1-w01)/wg1))) - np.sqrt(wg1-w01)*np.log(1-((compl1)/w01)**2)/(2*(compl1)**2))
        cp2 = -A2 * np.log(1-(compl2/w02)**2)/(2*(compl2)**2)
        return drude, cp1, cp2

    def _eps(self, omega, x_au):
        drude, cp1, cp2 = self._drude(omega, x_au)
        kk =  drude + cp1 + cp2
        eps1 = kk.real
        eps2 = kk.imag
        return eps1, eps2

    def get_n(self, wls):
        omega = h*c/wls
        eps1, eps2 = self._eps(omega, self.x_Au)
        n = (eps1 + 1j * eps2)**0.5
        return n.real

    def get_k(self, wls):
        omega = h*c/wls
        eps1, eps2 = self._eps(omega, self.x_Au)
        n = (eps1 + 1j * eps2)**0.5
        return n.imag

if __name__== '__main__':
    print('Alloy material test')
    from mstm_studio.alloy_AuAg import AlloyAuAg
    mat = AlloyAuAg(0.5)
    print(mat)
    # mat.plot()
    print('n, k = ', mat.get_n(800), mat.get_k(800))
    input('Press enter')
    mat.plot()
