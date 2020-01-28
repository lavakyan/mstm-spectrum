import setuptools

#with open("README.md", "r") as fh:
#    long_description = fh.read()


setuptools.setup(
    name="mst_studio-lavakyan",
    version="0.1",
    author="Leon Avakyan",
    author_email="laavakyan@sfedu.ru",
    description="Mie theory and T-matrix calculations GUI",
    #long_description=long_description,
    #long_description_content_type="text/markdown",
    long_description=long_description=
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
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
