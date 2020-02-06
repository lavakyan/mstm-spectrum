from mstm_studio.mstm_spectrum import Material
import numpy as np

wls = np.linspace(300, 800, 51)  # spectral region
omega = 1240. / wls              # freq. domian

omega_p = 9.        # plasma frequency, eV
gamma   = 0.03      # damping, eV
# Drude's dielectric function:
epsilon = 1 + omega_p**2 / (omega * (omega - 1j * gamma))

mat = Material('drude', eps=epsilon, wls=wls)
fig, axs = mat.plot()
fig.savefig('drude_eps.png', bbox_inches='tight')
