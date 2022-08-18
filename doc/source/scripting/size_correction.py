from mstm_studio.mstm_spectrum import Material
from mstm_studio.diel_size_correction import SizeCorrectedGold
from mstm_studio.contributions import MieSingleSphere
import numpy as np
import matplotlib.pyplot as plt

D = 3  # particle size
wls = np.linspace(400, 700, 201)  # spectral region

gold = Material('etaGold.txt')  # bulk dielectric function
gold_corr = SizeCorrectedGold('etaGold.txt')  # corrected, D-dependent

mie = MieSingleSphere(wavelengths=wls, name='MieSphere')

mie.set_material(material=gold, matrix=1.5)
plt.plot(wls, mie.calculate(values=[1, D]), '--', label='no size corr.')

mie.set_material(material=gold_corr, matrix=1.5)
plt.plot(wls, mie.calculate(values=[1, D]), label='size corrected')

plt.xlabel('Wavelength, nm')
plt.ylabel('Extinction efficiency')
plt.legend()
plt.savefig('size_correction.png')
plt.show()
