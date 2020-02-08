#!/usr/bin/python
# -*- coding: utf-8 -*-
from mstm_studio.alloy_AuAg import AlloyAuAg
from mstm_studio.contributions import LinearBackground, MieLognormSpheresCached
from mstm_studio.fit_spheres_optic import Fitter

fitter = Fitter(exp_filename='experiment.dat')        # load experiment from tabbed file
fitter.set_extra_contributions(
    [LinearBackground(fitter.wls),                    # wavelengths from experiment
     MieLognormSpheresCached(fitter.wls, 'LN Mie')],  # cached version for faster fittting
    [0.02, 0.0001, 0.1, 2.0, 0.4])                    # initial values for a, b, C, mu, sigma
fitter.extra_contributions[1].set_material(AlloyAuAg(1.), 1.5)  # gold particles in glass
fitter.set_spheres(None)  # no spheres - no slow MSTM runs

# run fit (takes ~20 seconds on 2GHz CPU)
fitter.run()

fitter.report_result()

# plot results
import matplotlib.pyplot as plt
plt.plot(fitter.wls, fitter.exp,  'ro', label='exp.')
plt.plot(fitter.wls, fitter.calc, 'b-', label='fit')
plt.xlabel('Wavelength, nm')
plt.ylabel('Exctinction, a.u.')
plt.legend()
plt.savefig('fit_by_Mie.png', bbox_inches='tight')
