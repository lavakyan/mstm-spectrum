#!/usr/bin/python
# -*- coding: utf-8 -*-
from mstm_studio.alloy_AuAg import AlloyAuAg
from mstm_studio.mstm_spectrum import SPR, ExplicitSpheres
import matplotlib.pyplot as plt
import numpy as np


wls = np.linspace(300, 800, 100)    # define wavelengths
spr = SPR(wls)                      # create SPR object
spr.environment_material = 'glass'  # set matrix

spheres = ExplicitSpheres(2, [0,0,0,10,0,0,0,12], mat_filename=[AlloyAuAg(1.),AlloyAuAg(0.)])
spr.set_spheres(spheres)
spr.set_incident_field(fixed=False)
spr.simulate()

plt.plot(spr.wavelengths, spr.absorbtion)
plt.xlabel('Wavelegth, nm')
plt.ylabel('Absorbtion efficiency, a.u.')
plt.savefig('core-shell_mstm.png', bbox_inches='tight')
