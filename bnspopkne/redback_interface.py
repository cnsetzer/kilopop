import numpy as np
from astropy import units as u
from bnspopkne.kne import saee_bns_emgw_with_viewing_angle as saee_pop_view


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
    wavelengths = (frequencies * u.Hz).to(u.Angstrom, equivalencies=u.spectral()).value
    kne = saee_pop_view(
        min_wave=wavelengths.min() - 500.0,  # do better kcorrection but for now do this
        max_wave=wavelengths.max(),
        sim_gw=False,
        EOS="sfho",
        EOS_path="/cfs/data/chse9649/input_data/kne_data/eos_data/",
        kappa_grid_path="/cfs/data/chse9649/input_data/kne_data/opacity_data/Setzer2022_popkNe_opacities.csv",
        gp_hyperparameter_file="/cfs/data/chse9649/input_data/kne_data/opacity_data/Setzer2022_popkNe_GP_hyperparam.npy",
        # EOS_path="/Users/cnsetzer/Documents/Project1_kNe/kne_modeling/eos_data/",
        # kappa_grid_path="Setzer2022_popkNe_opacities.csv",
        save=False,
        consistency_check=False,
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
