[![PyPI version](https://badge.fury.io/py/mstm-studio.svg)](https://badge.fury.io/py/mstm-studio) [![Documentation Status](https://readthedocs.org/projects/mstm-studio/badge/?version=latest)](https://mstm-studio.readthedocs.io/en/latest/?badge=latest)

# mstm-spectrum
## About
Python wrapper for multiple sphere T-matrix (MSTM) code and Mie theory to calculate surface plasmon resonance (SPR) spectrum and *fit* it to experiment.

Based and inspired by MSТM GUI code <https://git.stim.ee.uh.edu/optics/mstm-gui.git> by David Mayerich.

Multi-Sphere T-Matrix code is proposed by Dr. Daniel Mackowski and Dr. Michael Mishchenko in:
"A multiple sphere T-matrix Fortran code for use on parallel computer clusters,"
*Journal of Quantitative Spectroscopy and Radiative Transfer* (2011) 112 2182-2192
<https://dx.doi.org/10.1016/j.jqsrt.2011.02.019>.

Please cite the above reference if using MSTM code.

### Features

1. Materials defined from constant expression, from tabulated file or from analytical formula of Rioux et al [[doi](http://doi.org/10.1002/adom.201300457)] for Au-Ag;
1. Simple functional contributions (linear, lorentzian, gaussian)
1. Mie theory contributions
1. MSTM calculations
1. Fitting of experimental data by any of the mentioned contributions, including combinations
1. Interactive graphical user interface
1. Flexible Python scripting

### Installation

* source distribution is available through PyPi repository.
Install system-wide by command `pip3 install mstm_studio`
or for current user by command `pip3 install mstm_studio --user`.

* MSTM binary should be compiled and specified by environmental variable `MSTM_BIN`.
MSTM source and binaries can be obtained on [MSTM website](http://eng.auburn.edu/users/dmckwski/scatcodes/).
Precompiled binary for Linux Debian x64 and Windows x32 can be found in [latest release](releases/latest).

### Dependencies

* **Python** (tested: Python3 under Debian10 Linux; Anaconda Python3 under Windows7)
* **NumPy** - numerical python library
* **SciPy** - scientific python library

Optional
* **MatPlotLib** - plotting with python
* **tkinter**, **PIL** - tk libraries and Python image - for GUI
* **ScatterPy** - non-spherical particles (spheroids)

### Contributors

Avakyan L.A. <laavakyan@sfedu.ru>
and students Skidanenko A.V. (PhD) and Yablinovski K.A. (MS).


## Usage

### Graphical user interface

Under Linux can be run by `python3 -m mstm_studio` command.

Under Windows the following shell script may be used
```
@ECHO OFF
PATH=C:\ProgramData\Anaconda3;C:\ProgramData\Anaconda3\Library\mingw-w64\bin;C:\ProgramData\Anaconda3\Library\usr\bin;C:\ProgramData\Anaconda3\Library\bin;C:\ProgramData\Anaconda3\Scripts;C:\ProgramData\Anaconda3\bin;C:\ProgramData\Anaconda3\condabin;%PATH%
set MSTM_BIN="C:\Users\L\Desktop\mstm_studio old\mstm-spectrum\mstm.exe"
python.exe -m mstm_studio
PAUSE
```
In this script the `PATH` variable is updated to ensure python binary is in it.
Next, `MSTM_BIN` environmental variable is set.

Without MSTM the code will work in part, i.e. resctricted to Mie and simple functional (Gauss, Lorentz) contributions.

![GUI screenshot image][screen_gui]

### Python scripting

Alternatively, the python scripting way may be used, which
gives full control over the calculations.

Example script with fitting of experimental data by spheres of material with dielectric function specified in file.

``` python
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
plt.plot(fitter.wls, fitter.exp, 'ro',
         fitter.wls, fitter.calc, 'b-')
plt.xlabel('Wavelength, nm')
plt.ylabel('Exctinction, a.u.')
plt.show()
```

![Screenshot image][screen]


## Citation

If you used this code in a scientific paper please cite original MSTM code reference and the following reference:

L. Avakyan, M. Heinz, A. Skidanenko, K. Yablunovskiy, J. Ihlemann, J. Meinertz, C. Patzig, M. Dubiel, L. Bugaev
*J. Phys.: Condens. Matter* (2018) 30 045901 [[doi](http://doi.org/10.1088/1361-648X/aa9fcc)]


[screen_gui]: example/screenshot-gui.jpg?raw=true "GUI screenshot"
[screen]: example/screenshot-example.png?raw=true "Screenshot of example run"
