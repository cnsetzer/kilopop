import numpy as np
from astropy import units as u
from bnspopkne.kne import saee_bns_emgw_with_viewing_angle as saee_pop_view
from bnspopkne.macronovae_wrapper import make_rosswog_seds as mw
from bnspopkne import equation_of_state as eos
from bnspopkne import mappings, population
from sfdmap import SFDMap as sfd
from extinction import fitzpatrick99 as F99
from extinction import apply
from astropy.constants import c as speed_of_light_ms
from astropy.cosmology import Planck18 as cosmo
from redback.utils import calc_kcorrected_properties

# Set global module constants.
speed_of_light_kms = speed_of_light_ms.to("km/s").value  # Convert m/s to km/s
# Initialze sfdmap for dust corrections
sfdmap = sfd()


def redback_S22_BNS_popkNe(time, redshift, **kwargs):
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
    frequencies = kwargs["frequency"]
    kwargs["z"] = redshift
    freq_k, t_k = calc_kcorrected_properties(frequencies, redshift, time)
    wavelengths = (frequencies * u.Hz).to(u.Angstrom, equivalencies=u.spectral()).value
    kne = saee_pop_view(
        min_wave=wavelengths.min() - 500.0,  # do better kcorrection but for now do this
        max_wave=wavelengths.max(),
        **kwargs,
    )
    flux_density_rb_form = []
    for ti, wvi in list(zip(time, wavelengths)):
        flux_density_rb_form.append((kne.extincted_model.flux(ti, wvi)))
    # flux_density = (
    #     kne.extincted_model.flux(time_u, wave_flip)
    #     * u.erg
    #     / (u.cm ** 2 * u.s * u.Angstrom)
    # )
    # redback_flux = np.flip(flux_density, axis=1).reshape(-1)
    redback_flux = (
        np.asarray(flux_density_rb_form) * u.erg / (u.cm ** 2 * u.s * u.Angstrom)
    )
    if kwargs["output_format"] == "flux_density":
        return redback_flux.to(
            u.mJy, equivalencies=u.spectral_density(wavelengths * u.Angstrom)
        ).value
    elif kwargs["output_format"] == "magnitude":
        return redback_flux.to(u.ABmag).value


def redback_S22_BNS_popkNe_streamlined(
    time, redshift, m1, m2, theta_obs, disk_eff, peculiar_velocity, **kwargs
):
    """
    Streamlines the minimal necessary componenets needed to simulate the
    Setzer2022 et al. greybody population kilonova.

    :param time: observer frame time in days (observer-frame phase of the transient).
    :param redshift: redshift
    :param m1: The source-frame mass of the primary neutron star.
    :param m2: The source-frame mass of the secondary neutron star.
    :param theta_obs: The observer viewing angle of the kilonva.
    :param disk_eff: The disk unbinding efficiency for the secular ejecta.
    :param kwargs: output_format
                    frequency
    :return flux_density: can be in units of mJy or AB mag.
    """
    # m1 = kwargs["m1"]
    # m2 = kwargs["m2"]
    # theta_obs = kwargs["theta_obs"]
    # disk_eff = kwargs["disk_eff"]
    # peculiar_velocity = kwargs["peculiar_velocity"]

    frequencies = kwargs["frequency"]
    ra = kwargs["ra"]
    dec = kwargs["dec"]
    tov_mass = kwargs["tov_mass"]
    EOS_mass_to_rad = kwargs["EOS_mass_to_rad"]
    EOS_name = kwargs["EOS_name"]
    grey_opacity_interp = kwargs["grey_opacity_interp"]
    opacity_data = kwargs["opacity_data"]

    obs_z = (1 + redshift) * (
        np.sqrt(
            (1 + (peculiar_velocity / speed_of_light_kms))
            / ((1 - (peculiar_velocity / speed_of_light_kms)))
        )
    ) - 1.0
    freq_k, t_k = calc_kcorrected_properties(frequencies, obs_z, time)
    transient_duration = t_k.max()
    wavelengths_k = (freq_k * u.Hz).to(u.Angstrom, equivalencies=u.spectral()).value
    wavelengths = (frequencies * u.Hz).to(u.Angstrom, equivalencies=u.spectral()).value

    c1 = eos.compute_compactnesses_from_EOS(m1, EOS_mass_to_rad)
    c2 = eos.compute_compactnesses_from_EOS(m2, EOS_mass_to_rad)
    Ye = mappings.compute_ye_at_viewing_angle(theta_obs, EOS_name)
    m_ej_dyn, v_ej = mappings.map_to_dynamical_ejecta(m1, c1, m2, c2,)
    m_sec, m_tot, disk_eff = mappings.map_to_secular_ejecta(
        m1, c1, m2, c2, m_ej_dyn, tov_mass, disk_effs=disk_eff,
    )
    k_grey = mappings.map_kne_to_grey_opacity_via_gaussian_process(
        m_tot, v_ej, Ye, grey_opacity_interp, opacity_data,
    )

    KNE_parameters = []
    KNE_parameters.append(600.0 / 86400.0)  # starting time [days] 10 minutes
    KNE_parameters.append(transient_duration)  # ending time
    KNE_parameters.append(m_tot)  # total ejecta mass
    KNE_parameters.append(v_ej)  # median ejecta velocity
    KNE_parameters.append(1.3)  # nuclear heating rate exponent
    KNE_parameters.append(0.25)  # thermalization factor
    KNE_parameters.append(1.0)  # DZ enhancement
    KNE_parameters.append(k_grey)  # The grey opacity
    KNE_parameters.append(150.0)  # Initial temperature [K]
    KNE_parameters.append(Ye)  # Electron fraction
    # Not reading heating rates from file so feed fortran dummy
    # variables
    KNE_parameters.append(True)  # Flag to use numerical fit nuclear heating rates
    KNE_parameters.append(False)  # Read heating rates variable
    KNE_parameters.append("dummy string")  # Heating rates file
    _, _, flux = mw(KNE_parameters, times=t_k, wavelengths=wavelengths_k)
    dist_mpc = cosmo.luminosity_distance(redshift).value
    dist_pc = dist_mpc * 1000000.0  # convert Mpc to pc

    # Note that it is necessary to scale the amplitude relative to the 10pc
    # (i.e. 10^2 in the following eqn.) placement of the SED currently
    amp = np.power(10.0 / dist_pc, 2)
    flux *= amp
    # Current working around for issue with amplitude...

    r_v = 3.1
    unreddend_fluxes = flux
    reddend_fluxes = np.empty_like(unreddend_fluxes)
    uncorr_ebv = sfdmap.ebv(ra, dec, frame="icrs", unit="radian")
    for i, _ in enumerate(time):
        reddend_fluxes[i] = apply(
            F99(np.array([wavelengths[i]]), r_v * uncorr_ebv, r_v=r_v),
            unreddend_fluxes[i],
        )

    flux_density_rb_form = reddend_fluxes

    redback_flux = (
        np.asarray(flux_density_rb_form) * u.erg / (u.cm ** 2 * u.s * u.Angstrom)
    )

    if kwargs["output_format"] == "flux_density":
        output_flux = redback_flux.to(
            u.mJy, equivalencies=u.spectral_density(wavelengths * u.Angstrom)
        ).value
        return output_flux
    elif kwargs["output_format"] == "magnitude":
        return redback_flux.to(u.ABmag).value
