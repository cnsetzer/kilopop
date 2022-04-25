"""Functions to deal with the equation of state dependent properties of the model."""

import numpy as np
from scipy.interpolate import interp1d
from pandas import read_csv


def calculate_threshold_mass(tov_mass, EOS_mass_to_rad):
    """
    Function to calculate the prompt collapse threshold mass, given the
    maximum TOV mass and the radius at 1.6 solar mass for the chosen EOS.

    Returns:
    --------
        M_thr: float
            Threshold mass in solar masses.
    """
    a = 2.38
    b = 3.606
    M_tov = tov_mass
    R_16 = EOS_mass_to_rad(1.6)
    M_thr = (a - b * (M_tov / R_16)) * M_tov
    return M_thr


def get_EOS_table(EOS_path=None, EOS="sfho"):
    """
    Wrapper function to read the given Equation of State mass vs. radius
    diagram as pre-computed with a TOV-solver for use with this program.
    The format of the mass vs. radius data should be (radius, grav_mass,
    baryonic_mass,...)

    Parameters:
    -----------
        EOS_path: str
            The location of the file containing the mass vs. radius data.
        EOS: str
            The name of the equation of state.
        crust: boolean
            Flag indicating if the provided equation state contains a
            outer, stiff, neutron star crust.

    Returns:
    --------
        EOS_data: pandas.dataframe
            Dataframe containing the relationship between mass and radius
            for the given EOS.
    """
    EOS_data = read_csv(EOS_path + "mr_{}_full_right.csv".format(EOS))
    return EOS_data


def get_radius_from_EOS(EOS_table):
    """
    Wrapper function to create the interpolation function from the provided
    EOS data to evaluate the mass vs. radius relation for arbitray mass.

    Returns:
    --------
        f: scipy.interpolate instance
            Function which interpolates mass vs. radius for given EOS.
    """
    mass = EOS_table["grav_mass"]
    radius = EOS_table["radius"]
    f = interp1d(mass, radius)
    return f


def get_bary_mass_from_EOS(EOS_table):
    """
    Wrapper function to create the interpolation function from the provided
    EOS data to evaluate the mass vs. radius relation for arbitray mass.

    Returns:
    --------
        f: scipy.interpolate instance
            Function which interpolates mass vs. radius for given EOS.
    """
    mass = EOS_table["grav_mass"]
    bary_mass = EOS_table["bary_mass"]
    f = interp1d(mass, bary_mass)
    return f


def get_max_EOS_mass(EOS_table):
    """
    Wrapper function to find the Max TOV mass of the given EOS.

    Returns:
    --------
        tov_mass: float
            The maximum TOV mass.
    """
    tov_mass = max(EOS_table["grav_mass"].values)
    return tov_mass


def compute_compactnesses_from_EOS(mass, EOS_mass_to_rad):
    """
    Using the equation for compactness from neutron star mass and radius
    compute the compactnesses for both of the component stars in the system.

    Returns (Implicitly):
    ---------------------
        self.param3: float
            The stellar compactness of the first neutron star.
        self.param4: float
            The stellar compactness of the second neturon star.
    """
    G = 13.271317987 * np.power(10, 10)  # units km, M_sol^-1, (km/s)^2
    c = 299792.458  # units km/s
    R1 = EOS_mass_to_rad(mass)  # units km
    c1 = (G * mass) / ((c ** 2) * R1)
    return c1
