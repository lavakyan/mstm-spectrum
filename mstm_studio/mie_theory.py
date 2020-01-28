"""
  code from nanophotonics/npmie repository
  url: https://github.com/nanophotonics/npmie
  code by Alan Sanders

  modified to use in MSTM-studio by L.Aavakyan
"""

import numpy as np
from mstm_studio.mstm_spectrum import Material

try:
    from scipy.special import sph_jnyn
except:
    from scipy.special import spherical_yn, spherical_jn
    def sph_jnyn(maxn, z):
        jn  = []
        djn = []
        yn  = []
        dyn = []
        for n in range(0, maxn+1):
            jn.append (spherical_jn(n, z))
            djn.append(spherical_jn(n, z, derivative=True))
            yn.append (spherical_yn(n, z))
            dyn.append(spherical_yn(n, z, derivative=True))
        return np.array(jn), np.array(djn), np.array(yn), np.array(dyn)


def sph_hn(n, x):
    # calculate spherical hankel, h(n,x) = j(n,x) + iy(n,x) #
    jn, djn, yn, dyn = sph_jnyn(n, x)
    hn = jn + 1j * yn
    dhn = djn + 1j * dyn
    return hn, dhn


def calculate_mie_coefficients(n_max, x, m):
    """
    Calculates the Mie coefficients.

    :rtype : object
    :param n_max:
    :param x: size parameter
    :param m:
    """

    # calculate spherical bessels #
    jn, djn, yn, dyn = sph_jnyn(n_max, x)      # j(n, x), y(n, x)
    jm, djm, ym, dym = sph_jnyn(n_max, m * x)  # j(n, mx), y(n, mx)
    # calculate spherical hankel #
    hn, dhn = sph_hn(n_max, x)                 # h(n, x)
    # calculate riccati bessel functions #
    dpsi_n = [x * jn[n-1] - n * jn[n] for n in range(0, len(jn))]
    dpsi_m = [m * x * jm[n-1] - n * jm[n] for n in range(0, len(jm))]
    dzeta_n = [x * hn[n-1] - n * hn[n] for n in range(0, len(hn))]

    a_n = (m**2 * jm * dpsi_n - jn * dpsi_m) / (m**2 * jm * dzeta_n - hn * dpsi_m)
    b_n = (jm * dpsi_n - jn * dpsi_m) / (jm * dzeta_n - hn * dpsi_m)
    return a_n, b_n


def calculate_mie_efficiencies(r, wavelength, n_sph, n_med):
    """
    Calculates the mie efficiencies (q_scat, q_abs, q_ext, q_bscat)
    for a sphere in a dielectric medium at a given wavelength.

    :rtype : object
    :param r: radius of the sphere
    :param wavelength: wavelength of illumination
    :param n_sph: complex refractive index of the sphere
    :param n_med: real refractive index of the dielectric medium
    :return:
    """

    # calculate size parameter #
    x = n_med * (2 * np.pi / wavelength) * r    # x = n_med * kr, size parameter
    m = n_sph / n_med
    # n_max = int(np.ceil(x.real)+1)      # number of terms in series expansion
    n_max = int(x + 4 * x**(1.0 / 3.0) + 2)  # number of terms in series expansion

    q_scat = 0
    q_bscat = 0
    q_ext = 0
    q_abs = 0
    a_n, b_n = calculate_mie_coefficients(n_max, x, m)
    a = 0
    b = 0
    for n in range(1, n_max):
        a += a_n[n]
        b += b_n[n]
        q_scat += (2 * n + 1) * (abs(a_n[n])**2 + abs(b_n[n])**2)
        q_bscat += (2 * n + 1) * ((-1)**n) * (abs(a_n[n])**2 + abs(b_n[n])**2)
        q_ext += (2 * n + 1) * (a_n[n] + b_n[n]).real
    q_scat *= 2 / x**2
    q_bscat *= 2 / x**2
    q_ext *= 2 / x**2
    q_abs = q_ext - q_scat
    return q_scat, q_bscat, q_ext, q_abs


def calculate_mie_spectra(wavelengths, r, material, n_medium=1.):
    """
    Calculates the mie scattering and extinction efficiency of spherical
    nanoparticles with radius r and given material surrounded by a medium n_med
    for a set of given wavelengths.
    :rtype : object
    :param wavelengths: array of wavelengths to calculate spectra from
    :param r: radius of the sphere
    :param material: instance of Material class
    :param n_med: refractive index of the surrounding dielectric medium
    """

    mie_scattering = []
    mie_backscattering = []
    mie_extinction = []
    mie_absorption = []
    for wl in wavelengths:
        n_sph = material.get_n(wl) + 1j * material.get_k(wl)
        q_scat, q_bscat, q_ext, q_abs = calculate_mie_efficiencies(
            r, wl, n_sph, n_medium
        )
        mie_scattering.append(q_scat)
        mie_backscattering.append(q_bscat)
        mie_extinction.append(q_ext)
        mie_absorption.append(q_abs)
    return (np.array(mie_scattering), np.array(mie_backscattering),
            np.array(mie_extinction), np.array(mie_absorption))


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    #~ diameter_np = raw_input('Enter nanoparticle diameter (nm): ')
    #~ material = raw_input("Enter nanoparticle material: ")
    #~ medium = raw_input("Enter surrounding medium: ")
    diameter_np = material = medium = ''  # test
    if diameter_np == '':
        diameter_np = 140.
    else:
        diameter_np = float(diameter_np)
    if material == '':
        material = 'Au'
    if medium == '':
        medium = 1.
    else:
        medium = float(medium)
    mat_dict = {'Au': 'etaGold.txt', 'Ag': 'etaSilver.txt'}
    material_object = Material(3) # Material(mat_dict[material])
    wavelength = np.arange(300, 1000, 0.1)
    mie_scattering, mie_backscattering, mie_extinction, \
        mie_absorption = calculate_mie_spectra(
            wavelength, diameter_np / 2.0, material_object, medium
        )
    # save to file
    data = np.stack([wavelength, mie_scattering, mie_backscattering, \
        mie_extinction, mie_absorption])
    np.savetxt('MIE.dat', np.transpose(data), header='wl\tscatt\tbscatt\text\tabs')
    fig = plt.figure()
    # wavelength plots #
    ax = fig.add_subplot(411)
    ax.plot(wavelength, mie_scattering, 'r', label='scattering')
    ax.set_xticklabels(ax.get_xticklabels(), visible=False)
    ax.set_ylabel('scattering')
    ax = fig.add_subplot(412)
    ax.plot(wavelength, mie_backscattering, 'k', label='back-scattering')
    ax.set_xticklabels(ax.get_xticklabels(), visible=False)
    ax.set_ylabel('back-scattering')
    ax = fig.add_subplot(413)
    ax.plot(wavelength, mie_extinction, 'b', label='extinction')
    ax.set_xticklabels(ax.get_xticklabels(), visible=False)
    ax.set_ylabel('extinction')
    ax = fig.add_subplot(414)
    ax.plot(wavelength, mie_absorption, 'g', label='absorption')
    ax.set_ylabel('absorption')
    ax.set_xlabel('wavelength (nm)')
    plt.tight_layout()
    plt.show()
