#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
 Fit data in experiment.dat
 via Mie theory
"""

import mstm_studio as ms
from mstm_studio.fit_spheres_optic import Fitter
from mstm_studio.contributions import ConstantBackground, MieLognormSpheresCached
from mstm_studio.mstm_spectrum import Material


fitter = Fitter(exp_filename='experiment.dat')        # load experiment
fitter.set_extra_contributions(
    [ConstantBackground(fitter.wls),                  # specify wavelengths
     MieLognormSpheresCached(fitter.wls, 'LN Mie')],  # cached version for fittting
    [0.02, 0.1, 2.5, 0.2]                             # initial values for bkg, C, mu, sigma
    )
fitter.extra_contributions[1].set_material(Material('etaGold.txt'), 1.5)  # gold in glass
fitter.set_spheres(None)  # no spheres

# run fit (takes short time)
fitter.run()

fitter.report_result()

# plot results
import matplotlib.pyplot as plt
plt.plot(fitter.wls, fitter.exp, 'ro',
         fitter.wls, fitter.calc, 'b-')
plt.xlabel('Wavelength, nm')
plt.ylabel('Exctinction, a.u.')
plt.show()

