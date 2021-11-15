import numpy as np

from mstm_studio.mstm_spectrum import Material
from mstm_studio.contributions import MieSingleSphere


size_correction_gold = {'omp': 8.6, 'gbulk': 0.07, 'vF': 1.4, 'A': 0.3}


class SizeCorrectedMaterial(Material):

    def __init__(self, file_name, wls=None, nk=None, eps=None, sizecor=None):
        '''
        Create material with correction to the nanoparticle size
        (size induced limit to free path length).

        Parameters:

        file_name, wls, nk, eps : same meanining as for Material

        sizecor : dictionary for size correction parameters,
                  'omp' : plasma frequency in eV
                  'gbulk' : finite life-time broadening (bulk) in eV
                  'vF' : Fermi velocity (bulk) in nm*eV
                  'A' : unitless parameter

        size of the particle will be specified in self.D
        '''
        super().__init__(file_name, wls, nk, eps)
        self.D = None
        self.scprm = sizecor

    def _correction(self, wls):
        if self.D > 100:
            return 0
        # size is taken from self.D variable
        omp = self.scprm['omp']         # plasma freq. [eV]
        gbulk = self.scprm['gbulk']     # gamma bulk [eV]
        vF = self.scprm['vF'] / 1.6     # Fermi vel. [nm*eV]
        A = self.scprm['A']             # unitless constant
        omega = 1240 / wls  # nm -> eV
        gcorr = A * 2 * vF / self.D
        print('gamma correction: %f eV' % gcorr)
        return omp**2 / omega * ( 1 / (omega + 1j * gbulk) -
                1 / (omega + 1j * (gbulk + gcorr)))

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


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    wls = np.arange(300, 800, 1)

    matAu = Material('nk/etaGold.txt')

    matAu3nm = SizeCorrectedMaterial('nk/etaGold.txt', sizecor=size_correction_gold)
    if False:
        n = matAu.get_n(wls)
        k = matAu.get_k(wls)
        epsAu = (n + 1j * k)**2
        matAu3nm.set_size(3)
        n = matAu3nm.get_n(wls)
        k = matAu3nm.get_k(wls)
        epsAu3nm = (n + 1j * k)**2
        # ~ plt.plot(wls, matAu.get_n(wls), label='n')
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

    D = 2.17
    mie = MieSingleSphere(wavelengths=wls, name='ExtraContrib')
    mie.set_material(material=matAu, matrix=1.625)
    ext = mie.calculate(values=[1, D])

    mie = MieSingleSphere(wavelengths=wls, name='ExtraContrib')
    matAu3nm.D = D
    mie.set_material(material=matAu3nm, matrix=1.625)
    extcorr = mie.calculate(values=[1, D])


    plt.plot(wls, ext, label='%.0fnm, bulk' % D)
    plt.plot(wls, extcorr * np.sum(ext) / np.sum(extcorr), label='%.0fnm, corr.' % D)
    plt.xlabel('Wavelength, nm')
    plt.ylabel('Extinction')
    plt.legend()
    # ~ plt.savefig('Au%.0fnm_bulk_vs_sizecorr.png' % D)
    plt.show()


    print('See you!')

