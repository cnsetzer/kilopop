.. kilopop documentation master file, created by
   sphinx-quickstart on Mon Jan 23 14:13:46 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to kilopop's documentation!
================================================

:code:`kilopop` is a package to produce binary neutron star kilonovae in the grey-body approximation.
It can also create populations of these objects useful for forecasting detection and testing observing scenarios.
Additionally, it uses an emulator for the grey-opacity of the material calibrated against a suite of numerical radiation
transport simulations with the code :code:`SuperNu`. For more details on the components of this model see the accompanying paper
`Modelling Populations of Kilonovae <https://ui.adsabs.harvard.edu/abs/2022arXiv220512286S/abstract>`_.

.. automodule:: kilopop
    :members:

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   tutorial/tutorial_notebook_1.ipynb
   acknowledging

.. currentmodule:: kilopop


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


API:
----

.. autosummary::
   :toctree: api
   :template: custom-module-template.rst
   :caption: API:
   :recursive:

    equation_of_state
    kilonovae
    macronovae_wrapper
    mappings
    population_priors
