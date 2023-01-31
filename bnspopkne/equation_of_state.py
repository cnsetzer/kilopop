"""Module for equation of state dependent properties of the model."""

from astropy.constants import c, G
from scipy.interpolate import interp1d
from pandas import read_csv


def get_EOS_table(EOS_path=None):
    """Read in the Equation of state table.

    Wrapper function to read the given Equation of State mass vs. radius
    diagram as pre-computed with a TOV-solver for use with this program.
    The format of the mass vs. radius data should be (radius, grav_mass,
    baryonic_mass,...)

    Parameters
    -----------
    EOS_path: str
        The location of the file containing the mass vs. radius data. The
        provided data should be in the form of a labeled column csv.

    Returns
    --------
    EOS_data: pandas.dataframe
        Dataframe containing the relationship between mass and radius
        for the given EOS.
    """
    EOS_data = read_csv(EOS_path)
    return EOS_data


def get_radius_interpolator_from_EOS(EOS_table):
    """Read the radius table and create interpolator function.

    Wrapper function to create the interpolation function from the provided
    EOS data to evaluate the mass vs. radius relation for arbitray mass.

    Parameters
    -----------
    EOS_table: pandas DataFrame
        Data table of the gravitational mass vs. radius curve.

    Returns
    --------
    radius_interpolator: scipy.interpolate instance
        Function which interpolates mass vs. radius for given EOS.
    """
    mass = EOS_table["grav_mass"]
    radius = EOS_table["radius"]
    radius_interpolator = interp1d(mass, radius)
    return radius_interpolator


def get_max_EOS_mass(EOS_table):
    """Obtain the maximum supported neutron star mass from the given EOS.

    Wrapper function to find the Max TOV mass of the given EOS.

    Parameters
    -----------
    EOS_table: pandas DataFrame
        Data table of the gravitational mass vs. radius curve.

    Returns
    --------
    tov_mass: float
        The maximum mass supported by the given EOS.
    """
    tov_mass = max(EOS_table["grav_mass"].values)
    return tov_mass


def compute_compactnesses_from_EOS(mass, EOS_mass_to_rad):
    """Calucate the neutron star stellar compactness.

    Using the equation for compactness from neutron star mass and radius
    compute the compactnesses for both of the component stars in the system.

    Parameters
    -----------
    mass: float
        The gravitational mass in solar masses.
    EOS_mass_to_rad: function
        Interpolator function from provided EOS mass-radius curve.

    Returns
    ---------------------
    compactness: float
        The stellar compactness of the neutron star.
    """
    radius = EOS_mass_to_rad(mass)  # units km
    compactness = ((G.to('km3 / (M_sun s2)').value * mass) /
                   (((c.to('km/s').value) ** 2) * radius))
    return compactness
