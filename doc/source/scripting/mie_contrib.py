#!/usr/bin/python
# -*- coding: utf-8 -*-
from mstm_studio.contributions import MieLognormSpheres
from mstm_studio.alloy_AuAg import AlloyAuAg
import numpy as np

mie = MieLognormSpheres(name='mie', wavelengths=np.linspace(300,800,51))
mie.set_material(AlloyAuAg(x_Au=1), 1.5)  # golden sphere in glass

values = [1, 1.5, 0.5]  # scale, mu, sigma
fig, _ = mie.plot(values)
fig.savefig('mie_contrib.png', bbox_inches='tight')

mie.MAX_DIAMETER_TO_PLOT = 20   # 100 nm is by default
fig, _ = mie.plot_distrib(values)
fig.savefig('mie_distrib.png', bbox_inches='tight')
