'''
  code from nanophotonics/npmie repository
  url: https://github.com/nanophotonics/npmie
  code by Alan Sanders

  modified to use in MSTM-studio,
  added near-field by L.Aavakyan
'''

import numpy as np

from scipy.special import spherical_yn, spherical_jn, assoc_legendre_p
                                                      # lpmv ?


def _sph_jnyn(maxn, z):
    '''
    Calculate spherical Bessel functions
    jn, yn and thier derivatives
    using scipy.
    n = 0 .. maxn (including maxn)

    Parameters:
        maxn: int
            the highest rank

    Returns: 4 2D arrays of maxm, len(z) shape
        jn, djn, yn, dyn
        d means derivative
        as it was returned by obsolte `sph_jnyn` function
    '''
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


def _sph_hn(maxn, z):
    '''
    Returns: 2 2D arrays of maxm, len(z) shape
    hn, dhn -- function and it's derivative
    '''
    # calculate spherical hankel, h(n,x) = j(n,x) + iy(n,x) #
    jn, djn, yn, dyn = _sph_jnyn(maxn, z)
    hn = jn + 1j * yn
    dhn = djn + 1j * dyn
    return hn, dhn


def _vector_spherical_harmonics(n, m, theta, phi, x):
    '''
    Calculate vector spherical harmonics
    according to
    https://en.wikipedia.org/wiki/Vector_spherical_harmonics#Alternative_definition

    Input:
        n: int
            degree, n > 0
        m: int
            order, 0 < m < n+1. In our case of k||z always is m=1
        x: float
            x = k*r - size parameter
        theta: float or array-like
            inclination/polar angle (0 at z)
        phi: float or array-like
            azimutal angle (0 at x)

    Calc M N even/odd (1) / (3)
    projected on r, theta, and phi orts

    '''
    theta = np.asarray(theta)
    phi = np.asarray(phi)
    sinθ = np.sin(theta)
    cosθ = np.cos(theta)
    cosφ = np.cos(phi)
    sinφ = np.sin(phi)
    sinmθ = np.sin(m*theta)
    cosmθ = np.cos(m*theta)
    cosmφ = np.cos(m*phi)
    sinmφ = np.sin(m*phi)

    jn = spherical_jn(n, x)
    yn = spherical_yn(n, x)
    hn = jn + 1j * yn
    djn = spherical_jn(n, x, derivative=True)
    dyn = spherical_yn(n, x, derivative=True)
    dhn = djn + 1j * dyn
    dxjn = jn + djn  # d(x*jn(x))/dx
    dxhn = hn + dhn  # d(x*hn(x))/dx
    # pnm = lpmv(m=m, v=n, x=cosθ)  # associated Legendre polynomial
    #           assoc_legendre_p(n, m, z, ...
    pnm, dpnm = assoc_legendre_p(n, m, cosθ, branch_cut=2, # 2 or 3
                                 norm=False, diff_n=1)
    dpnm = -dpnm * sinθ

    m1_even_theta = - m / sinθ * sinmφ * pnm * jn
    m3_even_theta = - m / sinθ * sinmφ * pnm * hn
    m1_even_phi = - cosmθ * dpnm * jn
    m3_even_phi = - cosmθ * dpnm * hn

    m1_odd_theta = m / sinθ * cosmθ * pnm * jn
    m3_odd_theta = m / sinθ * cosmθ * pnm * hn
    m1_odd_phi = -sinmθ * dpnm * jn
    m3_odd_phi = -sinmθ * dpnm * hn

    n1_even_r = cosmφ * n * (n + 1) * pnm * jn / x
    n3_even_r = cosmφ * n * (n + 1) * pnm * hn / x
    n1_even_theta = cosmφ * dpnm * dxjn / x
    n3_even_theta = cosmφ * dpnm * dxhn / x
    n1_even_phi = -m / sinθ * sinmφ * pnm * dxjn / x
    n3_even_phi = -m / sinθ * sinmφ * pnm * dxhn / x

    n1_odd_r = sinmφ * n * (n + 1) * pnm * jn / x
    n3_odd_r = sinmφ * n * (n + 1) * pnm * hn / x
    n1_odd_theta = sinmφ * dpnm * dxjn / x
    n3_odd_theta = sinmφ * dpnm * dxhn / x
    n1_odd_phi = m / sinθ * cosmφ * pnm * dxjn / x
    n3_odd_phi = m / sinθ * cosmφ * pnm * dxhn / x

    return (m1_even_theta, m1_even_phi,
            m1_odd_theta,  m1_odd_phi,
            m3_even_theta, m3_even_phi,
            m3_odd_theta,  m3_odd_phi,
            n1_even_r, n1_odd_r,
            n1_even_theta, n1_even_phi,
            n1_odd_theta,  n1_odd_phi,
            n3_even_r, n3_odd_r,
            n3_even_theta, n3_even_phi,
            n3_odd_theta,  n3_odd_phi)

def _estimate_nmax(x):
    '''
    Estimate number of terms in series expansion nmax

    x: float, size parameter
    '''
    # n_max = int(np.ceil(x.real)+1)      # number of terms in series expansion
    return int(x + 4 * x**(1. / 3.) + 2)

def calculate_mie_coefficients(n_max, x, m, need_cd=False):
    '''
    Calculates the Mie coefficients.

    Parameters:
        n_max: int
            maximum expansion order
        x: float
            size parameter
        m: complex
            relative refr. index of the material
            m = n_sph / n_media
        need_cd: bool (False)
            True if needed coeff. `c` adn `d`
    Returns:
        a_n, b_n: arrays of Mie coefficients
    '''
    n_max = _estimate_nmax(x)
    jn, djn, yn, dyn = _sph_jnyn(n_max, x)      # j(n, x), y(n, x)
    jm, djm, ym, dym = _sph_jnyn(n_max, m * x)  # j(n, mx), y(n, mx)
    hn, dhn = _sph_hn(n_max, x)                 # h(n, x)
    # Riccati-Bessel functions:
    dpsi_n = [x * jn[n-1] - n * jn[n] for n in range(0, len(jn))]
    dpsi_m = [m * x * jm[n-1] - n * jm[n] for n in range(0, len(jm))]
    dzeta_n = [x * hn[n-1] - n * hn[n] for n in range(0, len(hn))]

    a_denom = m**2 * jm * dzeta_n - hn * dpsi_m  # common
    b_denom = jm * dzeta_n - hn * dpsi_m         # denominators
    a_n = (m**2 * jm * dpsi_n - jn * dpsi_m) / a_denom
    b_n = (jm * dpsi_n - jn * dpsi_m) / b_denom
    if need_cd:
        c_n = (jn * dzeta_n - hn * dpsi_n) / b_denom
        d_n = m * (jn * dzeta_n - hn * dpsi_n) / a_denom
        return a_n, b_n, c_n, d_n
    else:
        return a_n, b_n


def calculate_mie_efficiencies(r, wavelength, n_sph, n_med):
    '''
    Calculates the mie efficiencies (q_scat, q_abs, q_ext, q_bscat)
    for a sphere in a dielectric medium at a given wavelength.

    Parameters:
        r: float
            radius of the sphere
        wavelength: float
            wavelength of illumination
        n_sph: complex
            complex refractive index of the sphere
        n_med: float
            real refractive index of the dielectric medium

    Returns:
        q_scat, q_bscat, q_ext, q_abs
    '''
    x = n_med * (2 * np.pi / wavelength) * r  # n_med*k*r, size param
    n_max = _estimate_nmax(x)
    m = n_sph / n_med  # relative refr. index

    a_n, b_n = calculate_mie_coefficients(n_max, x, m)
    a_n = a_n[1:]  # n sarts from 1
    b_n = b_n[1:]  # a0 ~ 0 and b0 ~ 0 up to machine precision

    n = np.arange(1, n_max+1)
    weights = 2 * n + 1
    abs_sq = np.abs(a_n)**2 + np.abs(b_n)**2
    q_scat = np.sum(weights * abs_sq)
    sign = np.where(n % 2 == 0, 1, -1)
    # backscattering according to:
    # https://miepython.readthedocs.io/en/v2.2.2/07_algorithm.html
    q_bscat = np.abs(np.sum(weights * sign * (a_n - b_n)))**2
    q_ext = np.sum(weights * (a_n + b_n).real)

    q_scat *= 2 / x**2
    q_bscat *= 1 / x**2
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
        n_sph = material.get_nk(wl)
        q_scat, q_bscat, q_ext, q_abs = calculate_mie_efficiencies(
            r, wl, n_sph, n_medium
        )
        mie_scattering.append(q_scat)
        mie_backscattering.append(q_bscat)
        mie_extinction.append(q_ext)
        mie_absorption.append(q_abs)
    return (np.array(mie_scattering), np.array(mie_backscattering),
            np.array(mie_extinction), np.array(mie_absorption))


def calculate_mie_field(wavelength, r, material, n_medium=1.,
                        plane='XZ', umin=-20, umax=20,
                        vmin=-20, vmax=20, step=1):
    '''
    '''
    coef = n_medium * (2 * np.pi / wavelength)    # x = n_med * kr, size parameter
    m = material.get_nk(wavelength) / n_medium
    n_max = _estimate_nmax(coef * r)

    # prepare 2D grid
    us = np.arange(umin, umax, step)
    vs = np.arange(vmin, vmax, step)
    U, V = np.meshgrid(us, vs, indexing='xy')  # or indexing='ij'
    # map to spherical coordinates
    R = (U**2 + V**2)**0.5
    if plane.upper() in ['XZ', 'ZX', 'YZ', 'ZY']:
        phi = 0 if 'X' in plane.upper() else np.pi / 2.
        print(f'phi={phi}')
        Theta = np.arctan2(U, V)
    # ~ else:  # XY plane
        # ~ theta = np.pi / 2.
        # ~ Phi = np.arctan2(V, U)
        # ~ Phi[Phi < 0] += 2*np.pi

        sinθ = np.sin(Theta)
        cosθ = np.cos(Theta)
        cosφ = np.cos(phi)
        sinφ = np.sin(phi)

        def sph_to_cart(A_r, A_theta, A_phi):
            x = A_r * sinθ * cosφ + A_theta * cosθ * cosφ - A_phi * sinφ
            y = A_r * sinθ * sinφ + A_theta * cosθ * sinφ + A_phi * cosφ
            z = A_r * cosθ - A_theta * sinθ
            return x, y, z

        E_inc_r = np.zeros_like(R, dtype=complex)
        E_inc_theta = np.zeros_like(R, dtype=complex)
        E_inc_phi = np.zeros_like(R, dtype=complex)
        for n in range(1, n_max+1):
            m1_even_theta, m1_even_phi, \
            m1_odd_theta,  m1_odd_phi, \
            m3_even_theta, m3_even_phi, \
            m3_odd_theta,  m3_odd_phi, \
            n1_even_r, n1_odd_r, \
            n1_even_theta, n1_even_phi, \
            n1_odd_theta,  n1_odd_phi, \
            n3_even_r, n3_odd_r, \
            n3_even_theta, n3_even_phi, \
            n3_odd_theta,  n3_odd_phi = _vector_spherical_harmonics(
                                            n=n,
                                            m=1,
                                            theta=Theta,
                                            phi=phi,
                                            x=coef*R)

            weight = 1j**n * (2 * n + 1) / n / (n + 1)

            E_inc_r += -1j * weight * n1_even_r
            E_inc_theta += weight * (m1_odd_theta - 1j * n1_even_theta)
            E_inc_phi += weight * (m1_odd_phi - 1j * n1_even_phi)

        E_inc_x, E_inc_y, E_inc_z = sph_to_cart(E_inc_r, E_inc_theta, E_inc_phi)
        return (us, vs,
            E_inc_x, E_inc_y, E_inc_z
            )


if __name__ == '__main__':
    from mstm_studio.mstm_spectrum import Material
    import matplotlib.pyplot as plt
    diameter_np = 10.
    medium = 1.0
    material_object = Material(-1.0-0.1j) #
    # ~ material_object = Material('nk/etaSilver.txt')
    if True:
        wavelength = np.arange(300, 800, 1.)
        mie_scattering, mie_backscattering, mie_extinction, \
            mie_absorption = calculate_mie_spectra(
                wavelength, diameter_np / 2.0, material_object, medium
            )
        # save to file
        # ~ data = np.stack([wavelength, mie_scattering, mie_backscattering, \
            # ~ mie_extinction, mie_absorption])
        # ~ np.savetxt('MIE.dat', np.transpose(data), header='wl\tscatt\tbscatt\text\tabs')
        # wavelength plots #
        fig, axs = plt.subplots(4,1, sharex=True)
        ax = axs[0]
        ax.plot(wavelength, mie_scattering, 'r', label='scattering')
        ax.set_ylabel('scattering')
        ax = axs[1]
        ax.plot(wavelength, mie_backscattering, 'k', label='back-scattering')
        ax.set_ylabel('back-scattering')
        ax = axs[2]
        ax.plot(wavelength, mie_extinction, 'b', label='extinction')
        ax.set_ylabel('extinction')
        ax = axs[3]
        ax.plot(wavelength, mie_absorption, 'g', label='absorption')
        ax.set_ylabel('absorption')
        ax.set_xlabel('wavelength (nm)')
        plt.tight_layout()
        plt.show()
    if True:
        u, v, E_inc_x, E_inc_y, E_inc_z = calculate_mie_field(
            wavelength=200,
            r=diameter_np/2.,
            material=Material(-1.0-0.1j),
            n_medium=medium,
            plane='XZ',
            umin=-400, umax=400, vmin=-400, vmax=400, step=10)

        plt.pcolormesh(u, v, np.real(E_inc_x), cmap='RdBu')
        # ~ plt.pcolormesh(u, v, np.imag(E_inc_x), cmap='RdBu')
        # ~ E2 = np.abs(E_inc_x)**2 + np.abs(E_inc_y)**2 + np.abs(E_inc_z)**2
        # ~ plt.pcolormesh(u, v, E2, cmap='RdBu')
        plt.colorbar()
        plt.xlabel('X, nm')
        plt.ylabel('Z, nm')
        # ~ plt.legend()
        plt.tight_layout()
        plt.show()
