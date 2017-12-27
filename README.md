# mstm-spectrum
## About
Python wrapper for multiple sphere T-matrix (MSTM) code to calculate surface plasmon resonance (SPR) spectrum and *fit* it to experiment.

Based and inspired by MSТM GUI code <https://git.stim.ee.uh.edu/optics/mstm-gui.git> by David Mayerich.

Multi-Sphere T-Matrix code is proposed by Dr. Daniel Mackowski and Dr. Michael Mishchenko in:
"A multiple sphere T-matrix Fortran code for use on parallel computer clusters,"
*Journal of Quantitative Spectroscopy and Radiative Transfer* (2011) 112 2182-2192
<https://dx.doi.org/10.1016/j.jqsrt.2011.02.019>.

Please cite the above reference if using MSTM code.

## Usage

### Graphgical user interface

The intuitive (hope, it is) graphgical user interface
version can be run executing script `mstm_studio_support.py` or `mstm_studio.py`.

![GUI screenshot image][screen_gui]

### Python scripting

Alternatively, the python scripting way may choosen.
This is less intuitive, but allows to fine tune more details.

The possible workflow is:

1. Place files in the same directory:
    1. experiment file,
    1. dielectric constants files (i.e. `etaGold.txt`, `etaSilver.txt`)
    1. binaries named as:
        * for Windows: `mstm.exe`
        * for Linux: `run_mstm.sh` script to run `mstm.x` (see example)

        Source code and binaries can be obtained on [MSTM website](http://eng.auburn.edu/users/dmckwski/scatcodes/)
        or using [direct download link](http://eng.auburn.edu/users/dmckwski/scatcodes/mstm%20v3.0.zip).
        (You may contact us if suffer compilation problems)
1. Edit `start_fit.py` file to suit your needs. This will probably include:
    1. set to the directory with the scripts (remove these lines if scripts are stored in current directory):

        ``` python
        import sys
        # set the path to mstm_spectrum scripts. Binary mstm files should be in current folder.
        sys.path.append('/home/leon/ltg_projects/fit-T-matrix/mstm-spectrum')
        ```
    1. let python know whiсh classes we are going to use

        ``` python
        from mstm_spectrum import ExplicitSpheres
        from alloy_AuAg import AlloyAuAg
        from fit_spheres_optic import Fitter, FixConstraint
        ```
    1. set the experiment file name, background contribution ('constant', 'linear' or 'lorentz') and surrounding material:

        ``` python
        fitter = Fitter('example/optic_sample19.dat')
        fitter.set_background('lorentz')
        fitter.set_matrix('glass')  # 'glass', 'water', 'air' or explicit value, i.e. 1.66+0.1j
        ```
    1. specify the initial configuration of spheres.
        In this example, we will start from 4 spheres put in the verteces of square in XY plane:

      ``` python
        spheres = ExplicitSpheres(3, [-10, 10, 0], [0, 0, 14], [0, 0, 0],
                                  [12, 12, 12], [AlloyAuAg(x_Au=0.0)]*3)
        fitter.set_spheres(spheres)
        ```
    1. don't vary position of one of the spheres to save some computational resourses:

        ``` python
        fitter.add_constraint(FixConstraint('x0'))
        fitter.add_constraint(FixConstraint('y0'))
        fitter.add_constraint(FixConstraint('z0'))
        ```
    1. show some info, run, and report results:

        ``` python
        fitter.report_freedom()
        raw_input('Press enter')

        fitter.run()

        fitter.report_result()
        raw_input('press enter')
        ```
    1. execution of this script will show some stats:

        ```
        Number of spheres:      3
        Background parameters:  3
        Degree of freedom
        internal:       4
        external:       9

        ```
        After pressing enter (required bt `raw_input` command),
        the fitting will begin.
        At each accepted minimizing step the plot will be updated:

        ![Screenshot image][screen]


## Citation

If you used this code in a scientific paper please cite original MSTM code reference and the following reference:

L. Avakyan, M. Heinz, A. Skidanenko, K. Yablunovskiy, J. Ihlemann, J. Meinertz, C. Patzig, M. Dubiel, L. Bugaev
*J. Phys.: Condens. Matter* (2018) <http://doi.org/10.1088/1361-648X/aa9fcc>


[screen_gui]: images/screenshot-gui.jpg?raw=true "GUI screenshot"
[screen]: images/screenshot-example.jpg?raw=true "Screenshot of example run"
