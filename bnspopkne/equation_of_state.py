"""Module for equation of state dependent properties of the model."""

from astropy.constants import c, G
from scipy.interpolate import interp1d
from pandas import read_csv


def get_EOS_table(EOS_path=None):
    """
    Wrapper function to read the given Equation of State mass vs. radius
    diagram as pre-computed with a TOV-solver for use with this program.
    The format of the mass vs. radius data should be (radius, grav_mass,
    baryonic_mass,...)

    Parameters:
    -----------
        EOS_path: str
            The location of the file containing the mass vs. radius data.

    Returns:
    --------
        EOS_data: pandas.dataframe
            Dataframe containing the relationship between mass and radius
            for the given EOS.
    """
    EOS_data = read_csv(EOS_path)
    return EOS_data


def get_radius_from_EOS(EOS_table):
    """
    Wrapper function to create the interpolation function from the provided
    EOS data to evaluate the mass vs. radius relation for arbitray mass.
    Parameters:
    -----------
        EOS_table: pandas DataFrame
            Data table of the gravitational mass vs. radius curve.

    Returns:
    --------
        radius_interpolator: scipy.interpolate instance
            Function which interpolates mass vs. radius for given EOS.
    """
    mass = EOS_table["grav_mass"]
    radius = EOS_table["radius"]
    radius_interpolator = interp1d(mass, radius)
    return radius_interpolator


def get_max_EOS_mass(EOS_table):
    """
    Wrapper function to find the Max TOV mass of the given EOS.

    Parameters:
    -----------
        EOS_table: pandas DataFrame
            Data table of the gravitational mass vs. radius curve.

    Returns:
    --------
        tov_mass: float
            The maximum mass supported by the given EOS.
    """
    tov_mass = max(EOS_table["grav_mass"].values)
    return tov_mass


def compute_compactnesses_from_EOS(mass, EOS_mass_to_rad):
    """
    Using the equation for compactness from neutron star mass and radius
    compute the compactnesses for both of the component stars in the system.

    Parameters:
    -----------
        mass: float
            The gravitational mass in solar masses.
        EOS_mass_to_rad: function
            Interpolator function from provided EOS mass-radius curve.

    Returns:
    ---------------------
        compactness: float
            The stellar compactness of the neutron star.
    """
    radius = EOS_mass_to_rad(mass)  # units km
    compactness = ((G.to('km3 / (M_sun s2)').value * mass) /
                   (((c.to.('km/s').value) ** 2) * radius))
    return compactness


def calculate_threshold_mass(tov_mass, EOS_mass_to_rad):
    """
    Function to calculate the prompt collapse threshold mass, given the
    maximum TOV mass and the radius at 1.6 solar mass for the chosen EOS.
    Based on Bauswein et al. 2013.

    Parameters:
    -----------
        tov_mass: float
            The maximum mass from solving the TOV equations for given EOS.
        EOS_mass_to_rad: function
            Interpolator function from provided EOS mass-radius curve.
    Returns:
    --------
        prompt_collpase_mass_threshold: float
            Threshold mass in solar masses for prompt blackhole collapse.
    """
    # fit coefficients
    a = 2.38
    b = 3.606
    radius_1_6_km = EOS_mass_to_rad(1.6)
    prompt_collpase_mass_threshold = ((a - b *
                                      (tov_mass / radius_1_6_km)) * tov_mass)
    return prompt_collpase_mass_threshold
