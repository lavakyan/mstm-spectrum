import setuptools

#with open("README.md", "r") as fh:
#    long_description = fh.read()

package_data={'mstm_studio' : ['images/*16.png', 'images/splash.png']}

setuptools.setup(
    name="mstm_studio-lavakyan",
    version="0.2",
    author="Leon Avakyan",
    author_email="laavakyan@sfedu.ru",
    description="Mie theory and T-matrix calculations GUI",
    #long_description=long_description,
    #long_description_content_type="text/markdown",
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
    """,
    url="https://github.com/lavakyan/mstm-spectrum",
    packages=setuptools.find_packages(),
    install_requires=[
          'numpy',
          'scipy',
      ],
    extras_require={'GUI': ['matplotlib', 'tkinter', 'pillow']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    license='GPLv3',
    python_requires='>=3.6',
    include_package_data=True,
    package_data=package_data,
    #data_files=[('bitmaps', ['mstm_studio/images/splash.png', 'mstm_studio/images/sq_plus_icon&16.png'])],
)
