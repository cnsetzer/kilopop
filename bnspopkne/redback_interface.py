from astropy import units as u
from kne import saee_bns_emgw_with_viewing_angle as saee_pop_view


def redback_S22_BNS_popkNe(time, wavelengths: np.ndarray, output_unit:astropy.units.UnitBase = u.erg / (u.cm ** 2 * u.s * u.Angstrom), **kwargs):
    """
    Wrapper for the model to interface with redback fitting software.

    Parameters:
    -----------
        time: float
            The time at which to retrieve the flux density for the given transient
            instance.
        mass1: float
            The source frame mass of the primary component neutron star in solar
            masses.
        mass2: float
            The source frame mass of the secondary, i.e., smaller, neutron star
            in solar masses.
        inclination: float
            The inclination, and hence observer viewing angle (modulo axisymmetry)
            of the binary merger in radians.
        kwargs: dict
            Dictionary of other optional model parameters to be passed to the
            transient generator.
    """
    kne = saee_pop_view(kwargs)
    flux_density = kne.model.flux(time, wavelengths)*u.erg / (u.cm ** 2 * u.s * u.Angstrom)
    return flux_density.to(output_unit).value
