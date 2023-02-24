"""Public module for drawing population parameters."""
import numpy as np


def draw_disk_unbinding_efficiency(output_shape=None):
    """Set the disk unbinding efficiency.

    The unbinding efficiency is drawn from a uniform distribution over the
    range of values encountered in the literature of long-term GRMHD simulations
    of the remnant accretion disk around a neutron star merger remnant. See
    the paper for references to Metzger et al. 2008; Siegel & Metzger 2018;
    Miller et al. 2019; Fernández et al. 2019.

    Parameters
    -----------
    output_shape: float or tuple
        The shape/size of the efficiency values to draw.

    Returns
    --------
    disk_unbinding_efficiency: float or ndarray
        The fraction of matter unbound from the accretion disk contributing
        to the total radiating ejecta mass of the kilonova.
    """
    disk_unbinding_efficiency = np.random.uniform(0.1, 0.4, size=output_shape)
    return disk_unbinding_efficiency


def draw_viewing_angle(output_shape=None):
    """Set the observer viewing-angle.

    Function draw the observer angle at which the kilonovae is observed. This
    is drawn assuming uniform distribution of binary orbital plane alignment in
    the Universe and making the equivalence of the observer angle and the polar
    angle of the kNe, due to assumed axisymmetry about the normal axis to the
    binary merger plane.

    Parameters
    -----------
    output_shape: float or tuple
        The shape/size of the array to be created of viewing_angles

    Returns
    --------
    viewing_angle: float or ndarray
        The observer viewing angle, i.e., inclination, with respect to the
        binary merger plane.
    """
    # Viewing angle compute as a uniform distribution over the half-sphere due
    # to planar symmetry of the viewing angle in our model.
    viewing_angle = np.arccos(1.0 - np.random.random_sample(size=output_shape))
    # in radians
    return viewing_angle


def draw_mass_from_EOS_bounds(maximum_mass_bound, mass_lower_bound=1.0, output_shape=None):
    """Sample a uniform distribution for the mass.

    Parameters
    -----------
    maximum_mass_bound: float
        The maximum mass from which to draw the distribution. Typically,
        this is the TOV mass associated with the EOS. [solar masses]
    mass_lower_bound: float
        The minimum mass of the distribution, default 1.0 [solar mass].
    output_shape: float or tuple
        The number and array shaping structure for the draws.

    Returns
    --------
    source_frame_grav_mass: float
        The gravitational mass of the neutron star in the source frame.
    """
    source_frame_grav_mass = np.random.uniform(low=mass_lower_bound,
                                               high=maximum_mass_bound,
                                               size=output_shape)
    return source_frame_grav_mass


def draw_masses_from_EOS_bounds_with_mass_ratio_cut(
    maximum_mass_bound, mass_lower_bound=1.0, mass_ratio_cut=(2.0/5.0), output_shape=None
):
    """
    Draw neutron star component mass in the source frame given constraints.

    Draw the neutron star component masses assuming a uniform prior over
    the range of masses, as used by LIGO's BNS compact binary waveform search,
    but with a maximum mass set by the chosen EOS.

    Parameters
    -----------
    maximum_mass_bound: float
        The maximum TOV mass for the given EOS.
    mass_lower_bound: float
        The lower bound on the allowable mass of the neutron star.
    mass_ratio_cut: float
        The mass ratio constraint to apply to the mass sampling.
    output_shape: float or tuple
        The number and array shaping structure for the draws of masses.

    Returns
    ---------------------
    mass1: float
        The source-frame gravitational mass of the first neutron star
        [solar masses].
    mass2: float
        The source-frame gravitational mass of the second neutron star
        [solar masses].
    """
    mass1 = draw_mass_from_EOS_bounds(maximum_mass_bound,
                                      mass_lower_bound=mass_lower_bound,
                                      output_shape=output_shape)
    mass2 = draw_mass_from_EOS_bounds(mass1,
                                      mass_lower_bound=np.clip(mass1*mass_ratio_cut, a_min=1.0, a_max=None),
                                      output_shape=output_shape)
    return mass1, mass2
