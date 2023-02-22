kilopop

![](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/cnsetzer/f1a27976965673422ac94bc1afb240d3/raw/covbadge.json)


This package produces kilonovae following the population model presented in Setzer et al. 2023.

More Documentation is given here: https://cnsetzer.github.io/Setzer2022_BNSpopkNe/

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
