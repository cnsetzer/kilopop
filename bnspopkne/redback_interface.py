import numpy as np
from astropy import units as u
from kne import saee_bns_emgw_with_viewing_angle as saee_pop_view


def redback_S22_BNS_popkNe(
    time,
    redshift,
    frequency: np.ndarray,
    **kwargs
):
    """
    Wrapper for the model to interface with redback fitting software.

    Parameters:
    -----------
        time: float
            The time at which to retrieve the flux density for the given transient
            instance.
        redshift: float
            The redshift of the transients
        frequency: np.ndarray
            Frequencies to evaluate.
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
    # convert shapes
    time_u = np.unique(time)
    freq_u = np.unique(frequency)
    wavelengths = (freq_u*u.Hz).to(u.Ang, equivalenceies=u.spectral()).value
    kne = saee_pop_view(kwargs)
    flux_density = (
        kne.model.flux(time_u, wavelengths) * u.erg / (u.cm ** 2 * u.s * u.Angstrom)
    )
    flux_density_io = flux_density.to(output_unit).value
    redback_flux = flux_density_io.reshape(-1)
    if kwargs['output_format'] == 'flux_density':
        return (flux_density* u.erg / (u.cm ** 2 * u.s * u.Angstrom)).to(u.mJy, equivalenceies=u.spectral_density(wavelengths)).value
    elif kwargs['output_format'] == 'magnitude':
        return flux_density.to(u.ABmag).value
