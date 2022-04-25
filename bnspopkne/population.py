"""Public module for drawing population intrinsic parameters."""
import numpy as np


def draw_viewing_angle(inclinations=None, out_shape=1):
    """Set the observer viewing-angle.

    Function draw the observer angle at which the kilonovae is observed. This is
    drawn assuming uniform distribution of binary orbital plane alignment in the
    Universe and making the equivalence of the observer angle and the polar angle
    of the kNe, due to assumed axisymmetry about the normal axis to the binary
    merger plane.

    Parameters:
    -----------
        inclinations: ndarray (optional)
            Array of exisiting inclinations from which to replace certain values
            if set to None. Useful for redrawing a population if certain generated
            kNe were rejected due to parameter consistency checks.

    Returns:
    --------
        inclinations: float or ndarray
            The observer viewing angle, i.e., inclination, with respect to the
            binary merger plane.
    """
    if inclinations is None:
        o_shape = out_shape
    else:
        ind = np.argwhere(np.isnan(inclinations))
        o_shape = inclinations[ind].shape

    theta_obs = np.arccos(2 * np.random.random_sample(size=o_shape) - 1)  # in radians
    if inclinations is None:
        inclinations = theta_obs
    else:
        inclinations[ind] = theta_obs
    return inclinations


def draw_mass_from_EOS_bounds(max_mass, m_low=1.0, mass=None, out_shape=1):
    """Sample a uniform distribution for the mass.

    Parameters:
    -----------
        max_mass: float
            The maximum mass from which to draw the distribution.
        m_low: float
            The minimum mass of the distribution, default 1.0 solar mass.
        mass: float or nd.array
            The
    """
    if mass is None:
        ind = np.array([1])
        o_shape = out_shape
    else:
        ind = np.argwhere(np.isnan(mass))
        o_shape = mass[ind].shape

    new_mass = np.random.uniform(low=m_low, high=max_mass, size=o_shape)

    if mass is None:
        mass = new_mass
    else:
        mass[ind] = new_mass
    return mass


def draw_masses_from_EOS_bounds_with_mass_ratio_cut(max_mass, m_low=1.0, mass_ratio_cut=2.0 / 3.0, mass1=None, mass2=None, out_shape=1):
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
        mass1: nd.array-like
            The set of already provided masses, in solar masses, or default None.
        mass2: nd.array-like
            Same as mass1.
        out_shape: float or shape array
            The number and array shaping structure for the draws of masses.

    Returns (Implicitly):
    ---------------------
        m1: float
            The source-frame gravitational mass of the first neutron star [solar masses].
        m2: float
            The source-frame gravitational mass of the second neutron star [solar masses].
    """
    m1 = draw_mass_from_EOS_bounds(max_mass, mass=mass1, out_shape=out_shape)
    m2 = draw_mass_from_EOS_bounds(max_mass, mass=mass2, out_shape=out_shape)

    if np.isscalar(m1) is True:
        if m1 < m2:
            m1_exch = m1
            m1 = m2
            m2 = m1_exch
        while m2 / m1 < mass_ratio_cut:
            m1 = draw_mass_from_EOS_bounds(max_mass, mass=mass1, out_shape=out_shape)
            m2 = draw_mass_from_EOS_bounds(max_mass, mass=mass2, out_shape=out_shape)
            if m1 < m2:
                m1_exch = m1
                m1 = m2
                m2 = m1_exch
    else:
        ind3 = np.argwhere(m2 > m1)
        # mass ratio cut
        mass_q = np.divide(m2, m1)
        ind4 = np.argwhere(mass_q < mass_ratio_cut)
        ridx_1 = np.union1d(ind3, ind4)

        # To obtain uniform sampling in the m1,m2 plane resample if m2 > m1
        while ridx_1.shape[0] > 0:
            int_shape = np.shape(m2[ridx_1])
            m1[ridx_1] = draw_mass_from_EOS_bounds(max_mass, out_shape=int_shape)
            m2[ridx_1] = draw_mass_from_EOS_bounds(max_mass, out_shape=int_shape)
            ind3 = np.argwhere(m2 > m1)
            # mass ratio cut
            mass_q = np.divide(m2, m1)
            ind4 = np.argwhere(mass_q < mass_ratio_cut)
            # combine indicies into single set to resample
            ridx_1 = np.union1d(ind3, ind4)

    return m1, m2
