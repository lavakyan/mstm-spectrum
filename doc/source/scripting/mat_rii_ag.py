import numpy as np
import matplotlib.pyplot as plt
from mstm_studio.rii_materials import RiiMaterial

wls = np.arange(300, 800, 1)

mat = RiiMaterial()  # search RII archive in default paths
mat.scan()           # fill db items
print(mat.rii_db_items['main']['Ag'])
# ['Babar', 'Choi', 'Ciesielski', 'Ciesielski-Ge', 'Ciesielski-Ni',
#  'Ferrera-298K', 'Ferrera-404K', 'Ferrera-501K', 'Ferrera-600K',
#  'Hagemann', 'Jiang', 'Johnson', 'McPeak', 'Rakic-BB', 'Rakic-LD',
#  'Stahrenberg', 'Werner', 'Werner-DFT', 'Windt', 'Wu', 'Yang']

mat.select('main', 'Ag', 'Hagemann')
plt.plot(wls, mat.get_n(wls), 'b--', label='Ag-Hagemann')

mat.select('main', 'Ag', 'Jiang')
plt.plot(wls, mat.get_n(wls), 'm-.', label='Ag-Jiang')

mat.select('main', 'Ag', 'Johnson')
plt.plot(wls, mat.get_n(wls), 'r-', label='Ag-Johnson')


plt.legend()
plt.xlabel('Wavelength, nm')
plt.ylabel('Real part of complex refraction index')
plt.savefig('mat_rii_ag.png', bbox_inches='tight')
