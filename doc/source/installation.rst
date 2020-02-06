


Installation
============


Source code
-----------


Source code of Python wrapper is available on GitHub <https://github.com/lavakyan/mstm-spectrum>. 
Stable version published on PyPi <https://pypi.org/project/mstm-studio/>.

The source code of MSTM is not included and should be obtained from <https://scattport.org/index.php/light-scattering-software/multiple-particle-scattering/468-mstm>. 
MSTM studio can be run without MSTM binary, but with restricted functionality.


Linux installation
------------------


Install from PyPi:

``pip install mstm_studio``


On systems without root access:

``pip install mstm_studio --user``


Running GUI with 

``python -m mstm_studio``


May be required to explicitely specify python version, i.e. use ``pip3`` and ``python3`` in above commands.

Binding with MSTM
^^^^^^^^^^^^^^^^^

MSTM-studio will search for ``mstm.x`` binary in ``~/bin`` directory.
 
This can be altered by setting of `MSTM_BIN` environment variable, i.e. in bash:

``export $MSTM_BIN=~/my_compiled_mstm/mstm_v3.bin``


.. Note::   MSTM can be compiled with gfortran as::
      
       gfortran  mpidefs-serial.f90 mstm-intrinsics-v3.0.f90 mstm-modules-v3.0.f90 mstm-main-v3.0.f90 -O2  -o mstm.x
   
   This is serial compilation, for parallel the file ``mpidefs-serial.f90`` should be replaced. Consult the MSTM website for details.


Windows installation
--------------------

The tested way is using Anaconda Python distribution <https://www.anaconda.com/>. 

1. Open "Anaconda Prompt". The new terminal window should pop up. 
2. Type in ``pip install mstm_studio``. This may take a while since the dependent code will be downloaded and installed.
3. Check GUI by typing ``python -m mstm_studio`` in Anaconda Prompt 
   or check python scripting by typing ``import mstm_studio`` in python console.

Binding with MSTM
^^^^^^^^^^^^^^^^^

4. Obtain the MSTM binary. Put it to some folder. 
5. Setup environmental variable ``MSTM_BIN`` to point on the binary. 
   The shell comannd ``SETX MSTM_BIN="path_to_your_mstm_bin"`` 
   will do the temporary setup, which is useful for making ``*.cmd`` scripts. 
   Permanent setup of environemnt variable should be done with graphical interface, see for example, 
   <https://docs.oracle.com/en/database/oracle/r-enterprise/1.5.1/oread/creating-and-modifying-environment-variables-on-windows.html>.

.. Note:: If you write \*.cmd script to run gui, don't forget to update ``PATH`` variable to point on the Python distribution. 
    The easist way is to type ``echo %PATH%`` in Anaconda Promt, and use the output in your script.
    Example of GUI running script is ::
    
        @ECHO OFF
        PATH=C:\ProgramData\Anaconda3;C:\ProgramData\Anaconda3\Library\mingw-w64\bin;C:\ProgramData\Anaconda3\Library\usr\bin;C:\ProgramData\Anaconda3\Library\bin;C:\ProgramData\Anaconda3\Scripts;C:\ProgramData\Anaconda3\bin;C:\ProgramData\Anaconda3\condabin;%PATH%
        set MSTM_BIN="C:\Users\L\Desktop\mstm_studio old\mstm-spectrum\mstm.exe"
        python.exe -m mstm_studio
        PAUSE
        
    The last command (``PAUSE``) is put to prevent console windows from closing after program is ended.




