#!/usr/bin/python
# -*- coding: utf-8 -*-
from mstm_studio.alloy_AuAg import AlloyAuAg
from mstm_studio.contributions import LinearBackground, MieLognormSpheresCached
from mstm_studio.mstm_spectrum import ExplicitSpheres, Profiler
from mstm_studio.fit_spheres_optic import Fitter, FixConstraint, ConcentricConstraint

fitter = Fitter(exp_filename='experiment.dat')        # load experiment from tabbed file
fitter.set_extra_contributions(
    [LinearBackground(fitter.wls)],     # wavelengths from experiment
    [0.02, 0.0001])                     # initial values for a, b

spheres = ExplicitSpheres(2, [0,0,0,10,0,0,0,12], mat_filename=[AlloyAuAg(1.),AlloyAuAg(0.)])
fitter.set_spheres(spheres) # core-shell Au@Ag particle
fitter.set_matrix(1.5)      # in glass

fitter.add_constraint(ConcentricConstraint(0, 1))  # 0 -> 1
fitter.add_constraint(FixConstraint('x00'))
fitter.add_constraint(FixConstraint('y00'))
fitter.add_constraint(FixConstraint('z00'))

# run fit (takes ~200 seconds on 2GHz CPU)
with Profiler():
    fitter.run()

fitter.report_result()

# plot results
import matplotlib.pyplot as plt
plt.plot(fitter.wls, fitter.exp,  'ro', label='exp.')
plt.plot(fitter.wls, fitter.calc, 'b-', label='fit')
plt.xlabel('Wavelength, nm')
plt.ylabel('Exctinction, a.u.')
plt.legend()
plt.savefig('fit_by_core-shell.png', bbox_inches='tight')


# ~ ChiSq:      0.002354
# ~ Optimal parameters
        # ~ a00:        1.284882        (Varied:True)
        # ~ a01:        1.958142        (Varied:True)
        # ~ ext00:      0.186312        (Varied:True)
        # ~ ext01:      -0.000247       (Varied:True)
        # ~ scale:      -0.063814       (Varied:True)
        # ~ x00:        0.000000        (Varied:False)
        # ~ x01:        0.000000        (Varied:False)
        # ~ y00:        0.000000        (Varied:False)
        # ~ y01:        0.000000        (Varied:False)
        # ~ z00:        0.000000        (Varied:False)
        # ~ z01:        0.000000        (Varied:False)

