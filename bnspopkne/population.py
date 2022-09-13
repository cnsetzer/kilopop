"""Public module for drawing population intrinsic parameters."""
import numpy as np
from bnspopkne import equation_of_state as eos
from bnspopkne.kne import Setzer2022_kilonova
from bnspopkne import mappings


class Setzer2022_population_parameter_distribution(object):
    """
    Class to construct population from Setzer et al. 2022.

    This class holds attributes which are vectors of each parameter of the
    population.

    There is minimal tuning available to the user to generate a new population,
    otherwise the population from the paper is imported.
    """
    def __init__(
        self,
        population_size=50000,
    ):
        """
        Parameters:
        -----------
            population_size: int
                The size of the population to be drawn from the population.
        """
        self.population_size = population_size
        self.number_of_parameters = 12
        # Instantiate parameter values as none to be filled in later
        # Set the parameter names
        self.param1_name = "mass1"
        self.param2_name = "mass2"
        self.param3_name = "compactness1"
        self.param4_name = "compactness2"
        self.param5_name = "viewing_angle"
        self.param6_name = "electron_fraction"
        self.param7_name = "dynamical_ejecta_mass"
        self.param8_name = "median_ejecta_velocity"
        self.param9_name = "grey_opacity"
        self.param10_name = "secular_ejecta_mass"
        self.param11_name = "total_ejecta_mass"
        self.param12_name = "disk_unbinding_efficiency"

        # initialize attributes
        for i in range(self.number_of_parameters):
            setattr(self, "param{}".format(i + 1), np.empty(self.population_size))

        for i in range(self.population_size):
            kilonova = Setzer2022_kilonova(only_draw_parameters=True)
            for k in range(self.number_of_parameters):
                setattr(self, f"param{k + 1}"[i], getattr(kilonova, f"param{k + 1}"))


def draw_disk_unbinding_efficiency(output_shape=1):
    """Set the disk unbinding efficiency.

    Parameters:
    -----------
        output_shape: float or tuple
            The shape/size of the efficiency values to draw.
    Returns:
    --------
        disk_unbinding_efficiency: float or ndarray
            The fraction of matter unbound from the accretion disk contributing
            to the total radiating ejecta mass of the kilonova.
    """
    disk_unbinding_efficiency = np.random.uniform(0.1, 0.4, size=out_shape)
    return disk_unbinding_efficiency


def draw_viewing_angle(output_shape=1):
    """Set the observer viewing-angle.

    Function draw the observer angle at which the kilonovae is observed. This
    is drawn assuming uniform distribution of binary orbital plane alignment in
    the Universe and making the equivalence of the observer angle and the polar
    angle of the kNe, due to assumed axisymmetry about the normal axis to the
    binary merger plane.

    Parameters:
    -----------
        output_shape: float or tuple
            The shape/size of the array to be created of viewing_angles
    Returns:
    --------
        viewing_angle: float or ndarray
            The observer viewing angle, i.e., inclination, with respect to the
            binary merger plane.
    """
    viewing_angle = np.arccos(1.0 - np.random.random_sample(size=output_shape))
    # in radians
    return viewing_angle


def draw_mass_from_EOS_bounds(max_mass, m_low=1.0, output_shape=1):
    """Sample a uniform distribution for the mass.

    Parameters:
    -----------
        max_mass: float
            The maximum mass from which to draw the distribution.
        m_low: float
            The minimum mass of the distribution, default 1.0 solar mass.
        output_shape: float or tuple
            The number and array shaping structure for the draws.
    """
    source_frame_grav_mass = np.random.uniform(low=m_low,
                                               high=max_mass,
                                               size=output_shape)
    return source_frame_grav_mass


def draw_masses_from_EOS_bounds_with_mass_ratio_cut(
    max_mass, m_low=1.0, mass_ratio_cut=2.0 / 3.0, output_shape=1
):
    """
    Draw neutron star component mass in the source frame give constraints.

    Draw the neutron star component masses assuming a uniform prior over
    the range of masses, as used by LIGO's BNS compact binary waveform search,
    but with a maximum mass set by the chosen EOS.

    Parameters:
    -----------
        max_mass: float
            The maximum TOV mass for the given EOS.
        m_low: float
            The lower bound on the allowable mass of the neutron star.
        mass_ratio_cut: float
            The mass ratio constraint to apply to the mass sampling.
        output_shape: float or tuple
            The number and array shaping structure for the draws of masses.

    Returns (Implicitly):
    ---------------------
        m1: float
            The source-frame gravitational mass of the first neutron star
            [solar masses].
        m2: float
            The source-frame gravitational mass of the second neutron star
            [solar masses].
    """
    m1 = draw_mass_from_EOS_bounds(max_mass,
                                   m_low=m_low,
                                   output_shape=output_shape)
    m_low = m1*mass_ratio_cut
    m2 = draw_mass_from_EOS_bounds(max_mass,
                                   m_low=m_low,
                                   output_shape=output_shape)
    return m1, m2
