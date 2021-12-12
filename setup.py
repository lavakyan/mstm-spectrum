import setuptools

#with open("README.md", "r") as fh:
#    long_description = fh.read()

long_description= \
    """The graphical user interface and python scripts
    for calculation of optical extinction spectra of
    spherical particles used to study plasmonic
    nanoparticles.
    The binding to Multiple Spheres T-matrix (MSTM) code
    allows to simulate 'plasmonic' intercation between
    sphericl particles.
    Package allows basic analysis and fitting of experimental
    optical extinction spectra.
    """

package = 'mstm_studio'
version = __import__(package).__version__

package_data = {'mstm_studio' : ['images/*.png', 'nk/eta*.txt']}

setuptools.setup(
    name=package,
    version=version,
    author='Leon Avakyan',
    author_email='laavakyan@sfedu.ru',
    description='Mie theory and T-matrix calculations GUI',
    long_description=long_description,
    #long_description_content_type="text/markdown",
    url='https://github.com/lavakyan/mstm-spectrum',
    packages=setuptools.find_packages(),
    install_requires=[
          'numpy',
          'scipy',
      ],
    extras_require={#'MSTM' : ['tempfile'] ,
        'GUI': ['matplotlib', 'tkinter', 'pillow'],
        'spheroid': ['scatterpy']},
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    license='GPLv3',
    python_requires='>=3.8',
    include_package_data=True,
    package_data=package_data,
)
