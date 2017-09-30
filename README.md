# mstm-spectrum
## About
Python wrapper for multiple sphere T-matrix (MSTM) code to calculate surface plasmon resonance (SPR) spectrum

Based and inspired by MSMT GUI code <https://git.stim.ee.uh.edu/optics/mstm-gui.git> by David Mayerich.

Multi-Sphere T-Matrix code is proposed by Dr. Daniel Mackowski and Dr. Michael Mishchenko in:
Mackowski and Mishchenko, "A multiple sphere T-matrix Fortran code for use on parallel computer clusters," Journal of Quantitative Spectroscopy and Radiative Transfer, 112: 2182-2192, 2011.
<https://dx.doi.org/10.1016/j.jqsrt.2011.02.019>

Please cite the above reference if using MSTM code.

## Usage

1. Place files in the same directory:
    1. experiment file,
    1. dielectric constants files (i.e. `etaGold.txt`, `etaSilver.txt`)
    1. binaries named as:
        * for Windows: `mstm.exe`
        * for Linux: `run_mstm.sh` (to support parralel run via mpi) and `mstm.x`

        Source code and binaries can be obtained on [MSTM website](http://eng.auburn.edu/users/dmckwski/scatcodes/) or using [direct download link](<http://eng.auburn.edu/users/dmckwski/scatcodes/mstm%20v3.0.zip>).
1. Edit `start_fit.py` file to suit your needs. This will probably include:
    1. set to the directory with the scripts (remove this lines if they are stored in current directory):

        ```python
        import sys
        # set the path to mstm_spectrum scripts. Binary mstm files should be in current folder.
        sys.path.append('/home/leon/ltg_projects/fit-T-matrix/mstm-spectrum')
        ```
    1. set the experiment file name:

        ```python
        data = read_ascii('optic_sample22.dat', True, 0) # read and sort by 0th column
        ```
    1. set fitting interval and bins density:

        ```python
        wavelengths, exp = rebin(300, 800, 51, data[0,:], data[1,:])  # min 300 nm, max 800 nm, 51 bins
        ```

    1. set matrix material ('glass', 'water' and 'air' keywords are recognized) or refraction index value (1.0, 1.66, ..)

        ```
        fit_spheres_optic.MATRIX_MATERIAL = 'Glass'  # 1.66
        ```
    1. The values initialied before while-loop

        ```
        A = 60   # 'box' size
        a = 10   # sphere radius
        d = 10   # 'gap' between spheres
        ```
        are used to build initial configuration.
    1. at run time the script will show the number of spheres and prompt for continue, as shown below:

        ![Screenshot image][screen]

        In this example the first ~20 steps are preparational.
        They followed by production step, where chi-square is reported and plot is updated.

## Citation


[screen]: screenshot-example.jpg?raw=true "Screenshot of example run"
