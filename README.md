Setzer2022_BNSpopkNe

This package produces kilonovae following the population model presented in Setzer et al. 2022.

For basic usage, see example given in the Tutorial folder.

In order to execute this code the user will need to download and install
lapack 3.8 or higher from https://netlib.org/lapack/#_software.

Once this is done, the user must add the paths to these libraries to their local $LD_LIBRARY_PATH$ environment variable.

Though currently the working installation procedure is:

"""
python setup.py sdist

python setup.py build

python setup.py install
"""
