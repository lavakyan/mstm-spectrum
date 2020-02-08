#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
 Fit data in experiment.dat
 via MSTM theory
"""

from mstm_studio.fit_spheres_optic import Fitter
from mstm_studio.contributions import ConstantBackground
from mstm_studio.mstm_spectrum import ExplicitSpheres


fitter = Fitter(exp_filename='experiment.dat')   # load experiment
fitter.set_matrix(1.5)  # glass environment
fitter.set_extra_contributions([
    ConstantBackground(wavelengths=fitter.wls)
    ])
# initial configuration: three golden spheres at (0,0,0), (25,0,0) and (0,25,0) with radii 10.
spheres = ExplicitSpheres(N=3, Xc=[0,25,0], Yc=[0,0,25], Zc=[0,0,0], a=[10,10,10], mat_filename='etaGold.txt')
fitter.set_spheres(spheres)

# run fit (takes about a hour)
fitter.run()

fitter.report_result()


# plot results
import matplotlib.pyplot as plt
plt.plot(fitter.wls, fitter.exp, 'ro', label='exp.')
plt.plot(fitter.wls, fitter.calc, 'b-', label='fit')
plt.xlabel('Wavelength, nm')
plt.ylabel('Exctinction, a.u.')
plt.legend()
plt.show()



#ChiSq:  0.000572
#Optimal parameters
#        a00:    12.024634       (Varied:True)
#        a01:    10.268477       (Varied:True)
#        a02:    10.388544       (Varied:True)
#        ext00:  0.005471        (Varied:True)
#        scale:  0.013125        (Varied:True)
#        x00:    1.147471        (Varied:True)
#        x01:    25.361571       (Varied:True)
#        x02:    -5.158585       (Varied:True)
#        y00:    0.350451        (Varied:True)
#        y01:    1.000000        (Varied:True)
#        y02:    25.569792       (Varied:True)
#        z00:    1.000000        (Varied:True)
#        z01:    1.000000        (Varied:True)
#        z02:    1.000000        (Varied:True)

