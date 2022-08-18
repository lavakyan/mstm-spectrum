import numpy as np
from mstm_studio.mstm_spectrum import Material


class SizeCorrectedMaterial(Material):

    def __init__(self, file_name, wls=None, nk=None, eps=None,
                 omega_p=10.0, gamma_b=0.1, v_Fermi=1.0, sc_C=0.0):
        '''
        Create material with correction to the finite crystal size.
        Only life-time limit (~1/gamma) is considered.
        This should be sufficient for the sizes above ~2 nm.
        The particles smaller than ~2 nm require more sofisticated
        modifications (band gap, etc.)

        Parameters:

        file_name, wls, nk, eps:
            same meanining as for Material

        Parameters for size correction:

        omega_p:
            plasma frequency (bulk) [eV]

        gamma_b:
            life-time broadening (bulk) [eV]

        v_Fermi:
            Fermi velocity (bulk) [nm/fs]

        sc_C:
            size-corr. adj. parameter [unitless]

        size of the particle is specified as `self.D`
        '''
        super().__init__(file_name, wls, nk, eps)

        self.omega_p = omega_p
        self.gamma_b = gamma_b
        self.v_Fermi = v_Fermi
        self.sc_C = sc_C
        if np.abs(sc_C) < 1e-3:
            print('WARNING: Size correction parameters are not set')

        self.D = 1000

    def get_gamma_corr(self):
        '''
        correction to the life time energy broadening
        (gamma) in the Drude low induced by the
        finite particle size
        '''
        vF = self.v_Fermi / 1.6  # Fermi vel. in [nm*eV]
        return self.sc_C * 2 * vF / self.D

    def _correction(self, wls):
        # size is taken from self.D variable
        if self.D > 100:  # don't consider very big particles
            return 0
        omega = 1240 / wls  # nm -> eV
        # correction to gamma due to mean free path limit by size
        gamma_corr = self.get_gamma_corr()
        return self.omega_p**2 / omega * \
            (1 / (omega + 1j * self.gamma_b) -
             1 / (omega + 1j * (self.gamma_b + gamma_corr)))

    def set_size(self, D):
        if D > 0:
            self.D = D
        else:
            self.D = None

    def get_n(self, wls):
        if self.D is None:
            return super().get_n(wls)
        n = super().get_n(wls)
        k = super().get_k(wls)
        eps = (n + 1j*k)**2
        eps += self._correction(wls)
        mod = np.absolute(eps)
        n = np.sqrt((mod + np.real(eps)) / 2)
        return n

    def get_k(self, wls):
        if self.D is None:
            return super().get_k(wls)
        n = super().get_n(wls)
        k = super().get_k(wls)
        eps = (n + 1j*k)**2
        eps += self._correction(wls)
        mod = np.absolute(eps)
        k = np.sqrt((mod - np.real(eps)) / 2)
        return k


class SizeCorrectedGold(SizeCorrectedMaterial):

    def __init__(self, file_name, wls=None, nk=None, eps=None):
        '''
        Size correction for gold dielectric function
        (mean free path is limited by particle size)
        according to

        A. Derkachova, K. Kolwas, I. Demchenko
        Plasmonics, 2016, 11, 941
        doi: <10.1007/s11468-015-0128-7>
        '''
        super().__init__(file_name, wls, nk, eps,
                         omega_p=8.6,
                         gamma_b=0.07,
                         v_Fermi=1.4,
                         sc_C=0.3)


class SizeCorrectedSilver(SizeCorrectedMaterial):

    def __init__(self, file_name, wls=None, nk=None, eps=None):
        '''
        Size correction for gold dielectric function
        (mean free path is limited by particle size)
        according to

        J.M.J. Santillán, F.A. Videla, M.B.F. van Raap, D. Muraca,
        L.B. Scaffardi, D.C. Schinca
        J. Phys. D: Appl. Phys., 2013, 46, 435301
        doi: <10.1088/0022-3727/46/43/435301>
        '''
        super().__init__(file_name, wls, nk, eps,
                         omega_p=8.9,
                         gamma_b=0.111,
                         v_Fermi=1.41,
                         sc_C=0.8)


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from mstm_studio.contributions import MieSingleSphere

    wls = np.arange(300, 800, 1)

    ### Gold ###
    matAu = Material('nk/etaGold.txt')
    matAu3nm = SizeCorrectedGold('nk/etaGold.txt')

    if False:
        n = matAu.get_n(wls)
        k = matAu.get_k(wls)
        epsAu = (n + 1j * k)**2
        matAu3nm.set_size(3)
        n = matAu3nm.get_n(wls)
        k = matAu3nm.get_k(wls)
        epsAu3nm = (n + 1j * k)**2
        # plt.plot(wls, matAu.get_n(wls), label='n')
        plt.plot(1240/wls, np.imag(epsAu), label='bulk')
        plt.plot(1240/wls, np.imag(epsAu3nm), label='3nm')
        plt.xlabel('E, eV')
        plt.ylabel('Im eps')
        plt.legend()
        plt.show()
        plt.plot(1240/wls, np.real(epsAu), label='bulk')
        plt.plot(1240/wls, np.real(epsAu3nm), label='3nm')
        plt.xlabel('E, eV')
        plt.legend()
        plt.show()

        plt.plot(wls, matAu.get_n(wls), label='bulk')
        plt.plot(wls, matAu3nm.get_n(wls), label='3nm')
        plt.xlabel('Wavelength, nm')
        plt.ylabel('n')
        plt.legend()
        plt.show()
        plt.plot(wls, matAu.get_k(wls), label='bulk')
        plt.plot(wls, matAu3nm.get_k(wls), label='3nm')
        plt.xlabel('Wavelength, nm')
        plt.ylabel('k')
        plt.show()

    D = 3
    mie = MieSingleSphere(wavelengths=wls, name='ExtraContrib')
    mie.set_material(material=matAu, matrix=1.5)
    ext = mie.calculate(values=[1, D])

    mie = MieSingleSphere(wavelengths=wls, name='ExtraContrib')
    matAu3nm.D = D
    mie.set_material(material=matAu3nm, matrix=1.5)
    extcorr = mie.calculate(values=[1, D])

    plt.plot(wls, ext, label='%.0fnm, bulk' % D)
    plt.plot(wls, extcorr,  # * np.sum(ext) / np.sum(extcorr),
             label='%.0fnm, corr.' % D)
    plt.xlabel('Wavelength, nm')
    plt.ylabel('Extinction')
    plt.legend()
    plt.savefig('Au%.0fnm_bulk_vs_sizecorr.png' % D)
    plt.show()

    ### Silver ###
    matAg = Material('nk/etaSilver.txt')
    matAg3nm = SizeCorrectedSilver('nk/etaSilver.txt')

    mie = MieSingleSphere(wavelengths=wls, name='ExtraContrib')
    mie.set_material(material=matAg, matrix=1.5)
    ext = mie.calculate(values=[1, D])

    mie = MieSingleSphere(wavelengths=wls, name='ExtraContrib')
    matAg3nm.D = D
    mie.set_material(material=matAg3nm, matrix=1.5)
    extcorr = mie.calculate(values=[1, D])

    plt.plot(wls, ext, label='%.0fnm, bulk' % D)
    plt.plot(wls, extcorr,  # * np.sum(ext) / np.sum(extcorr),
             label='%.0fnm, corr.' % D)
    plt.xlabel('Wavelength, nm')
    plt.ylabel('Extinction')
    plt.legend()
    plt.savefig('Ag%.0fnm_bulk_vs_sizecorr.png' % D)
    plt.show()

    print('See you!')
