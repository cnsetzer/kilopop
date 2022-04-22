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
