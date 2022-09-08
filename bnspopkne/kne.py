"""This is a sample docstring for module."""
import warnings
from copy import deepcopy
import numpy as np
from astropy.constants import c as speed_of_light_ms
import astropy.units as units
from astropy.time import Time
from sncosmo import TimeSeriesSource, Model
from astropy.cosmology import Planck18 as cosmos
from astropy.cosmology import z_at_value
from sfdmap import SFDMap as sfd
from extinction import fitzpatrick99 as F99
from extinction import apply
from .macronovae_wrapper import make_rosswog_seds as mw
from bnspopkne.inspiral import compact_binary_inspiral
from bnspopkne import equation_of_state as eos
from bnspopkne import mappings
from bnspopkne import population
from scipy.integrate import quadrature

warnings.filterwarnings("ignore", message="ERFA function")

# Set global module constants.
speed_of_light_kms = speed_of_light_ms.to("km/s").value  # Convert m/s to km/s
# Initialze sfdmap for dust corrections
sfdmap = sfd()


class em_transient(object):
    """
    Base class for transient instances.

    This groups class methods that should be common to all transient models.
    These include assign the transient a peculiar velocity, redshifting, etc.

    """

    def __init__(self, z=0.001, cosmo=cosmos):
        """Init of the em_transient class wrapper."""
        source = TimeSeriesSource(self.phase, self.wave, self.flux)
        self.model = deepcopy(Model(source=source))
        # Use deepcopy to make sure full class is saved as attribute of new class
        # after setting model empty phase, wave, flux attributes
        self.phase = None
        self.wave = None
        self.flux = None
        self.put_in_universe(z, cosmo)

    def put_in_universe(self, z, cosmo=cosmos, pec_vel=None, dl=None, r_v=3.1):
        """Place transient instance into the simulated Universe.

        It sets spacetime location, ra, dec, t0, redshift, and properties of the
        transient due to the environment such as peculiar velocity and dervied
        observed redshift.

        Input parameters:
        -----------------

        """
        if z and not dl:
            self.z = z
            self.dist_mpc = cosmo.luminosity_distance(z).value
        elif dl and not z:
            self.dist_mpc = dl
            self.z = z_at_value(cosmo.luminosity_distance, dl * units.Mpc).value
        elif z and dl:
            self.dist_mpc = dl
            self.z = z
        else:
            raise ValueError

        self.dist_pc = self.dist_mpc * 1000000.0  # convert Mpc to pc
        self.peculiar_velocity(pec_vel)
        self.redshift()
        self.tmax = self.t0 + self.model.maxtime()
        self.extinct_model(r_v=3.1)
        if self.save is True:
            self.save_info(cosmo)

    def redshift(self):
        """Redshift the source.

        Wrapper function to redshift the spectrum of the transient instance,
        and scale the flux according to the luminosity distance.

        Input Parameters:
        -----------------
            cosmo: Astropy.cosmology instance
                Class instance of a cosmology that groups the functions needed
                to compute cosmological quantities.

        Output:
        -------
            self: modified self class instance
        """
        self.model.set(z=self.obs_z)
        # Note that it is necessary to scale the amplitude relative to the 10pc
        # (i.e. 10^2 in the following eqn.) placement of the SED currently
        # lumdist = self.dist_mpc * 1e6  # in pc
        # amp = np.power(10.0 / lumdist, 2)
        # self.model.set(amplitude=amp)

        # Current working around for issue with amplitude...
        mapp = 5.0 * np.log10(self.dist_pc / 10.0) + self.model.source_peakmag(
            "lsstz", "ab", sampling=0.05
        )
        self.model.set_source_peakmag(m=mapp, band="lsstz", magsys="ab", sampling=0.05)

    def extinct_model(self, r_v=3.1):
        """Apply dust extinction to transient."""
        phases = np.linspace(self.model.mintime(), self.model.maxtime(), num=1000)
        waves = np.linspace(self.model.minwave(), self.model.maxwave(), num=2000)
        unreddend_fluxes = self.model.flux(phases, waves)
        reddend_fluxes = np.empty_like(unreddend_fluxes)
        uncorr_ebv = sfdmap.ebv(self.ra, self.dec, frame="icrs", unit="radian")
        for i, phase in enumerate(phases):
            reddend_fluxes[i, :] = apply(
                F99(waves, r_v * uncorr_ebv, r_v=r_v), unreddend_fluxes[i, :]
            )

        source = TimeSeriesSource(phases, waves, reddend_fluxes)
        self.extincted_model = deepcopy(Model(source=source))

    def peculiar_velocity(self, pec_vel=None):
        """Set peculiar velocity.

        Draw from Gaussian peculiar velocity distribution with width 300km/s
        Reference: Hui and Greene (2006) Also, use transient id to set seed for
        setting the peculiar velocity, this is for reproducibility.
        """
        state = np.random.get_state()
        np.random.seed(seed=self.id)
        if pec_vel is None:
            # Apply typical peculiar_velocity type correction
            self.peculiar_vel = np.random.normal(loc=0, scale=300)
        else:
            self.peculiar_vel = pec_vel
        np.random.set_state(state)
        self.obs_z = (1 + self.z) * (
            np.sqrt(
                (1 + (self.peculiar_vel / speed_of_light_kms))
                / ((1 - (self.peculiar_vel / speed_of_light_kms)))
            )
        ) - 1.0

    def save_info(self, cosmo):
        """
        Basic basic docstring.
        """
        lsst_bands = ["lsstu", "lsstg", "lsstr", "lssti", "lsstz", "lssty"]
        times = np.linspace(0.0, self.model.maxtime(), 1001)
        for band in lsst_bands:
            setattr(
                self,
                f"peak_{band}",
                self.extincted_model.source_peakmag(band, "ab", sampling=0.1),
            )
            setattr(
                self,
                f"peak_abs_{band}",
                self.model.source_peakabsmag(band, "ab", sampling=0.1, cosmo=cosmo),
            )
            peak_mag = getattr(self, f"peak_{band}")
            lc_mags = self.model.bandmag(band, "ab", time=times)

            one_mag_inds = np.nonzero(lc_mags <= peak_mag + 1)
            one_mag_total_time = times[one_mag_inds][-1] - times[one_mag_inds][0]
            peak_time_index = np.argmin(lc_mags)
            peak_time = times[peak_time_index]
            day1_ind = np.argmin(abs(times - peak_time - 1))
            day2_ind = np.argmin(abs(times - peak_time - 2))
            day3_ind = np.argmin(abs(times - peak_time - 3))
            day4_ind = np.argmin(abs(times - peak_time - 4))
            day5_ind = np.argmin(abs(times - peak_time - 5))
            full_ind = -1
            dmdt_1 = -(peak_mag - lc_mags[day1_ind]) / 1.0
            dmdt_2 = -(peak_mag - lc_mags[day2_ind]) / 2.0
            dmdt_3 = -(peak_mag - lc_mags[day3_ind]) / 3.0
            dmdt_4 = -(peak_mag - lc_mags[day4_ind]) / 4.0
            dmdt_5 = -(peak_mag - lc_mags[day5_ind]) / 5.0
            dmdt_full = -(peak_mag - lc_mags[full_ind]) / (times[full_ind] - peak_time)
            setattr(self, f"peak_delay_from_merger_{band}", peak_time)
            setattr(self, f"onemag_peak_duration_{band}", one_mag_total_time)
            setattr(self, f"oneday_postpeak_dmdt_{band}", dmdt_1)
            setattr(self, f"twoday_postpeak_dmdt_{band}", dmdt_2)
            setattr(self, f"threeday_postpeak_dmdt_{band}", dmdt_3)
            setattr(self, f"fourday_postpeak_dmdt_{band}", dmdt_4)
            setattr(self, f"fiveday_postpeak_dmdt_{band}", dmdt_5)
            setattr(self, f"full_postpeak_dmdt_{band}", dmdt_full)


class kilonova(em_transient, compact_binary_inspiral):
    """Base class for kilonova transients, groups relevant class methods and attributes."""

    def __init__(self, z, cosmo, sim_gw):
        """Init class.

        Wrapper class to handle different kilonva models.

        """
        self.type = "kne"
        em_transient.__init__(self, z, cosmo)
        if sim_gw is True:
            compact_binary_inspiral.__init__(self)


class saee_bns_emgw_with_viewing_angle(kilonova):
    """
    Top-level class for kilonovae transients based on Rosswog, et. al 2017
    semi-analytic model for kilonovae spectral energy distributions.

    Parameters:
    -----------
        m1: float
            The gravitational mass of the first neutron star.
        m2: float
            The gravitational mass of the second neutron star.
        c1: float
            The stellar compactness of the first neutron star.
        c2: float
            The stellar compactness of the second neutron star.
        theta_obs: float
            The viewing angle of the observer with respect to the merger.
        Y_e: float
            The electron fraction along the line of sight of the observer.
        m_ej_dyn: float
            The dynamic ejecta mass of the expanding kilonova material.
        vej: float
            The mean ejecta velocity of the expanding kilonova material.
        kappa: float
            The grey opacity of the material along the line of sight of the
            observer.
        m_ej_sec: float
            The secular ejecta mass of the expanding kilonova material.
        m_ej_total: float
            The total ejecta mass of the expanding kilonova material.
        EOS: str
            The name of the equation of state of the neutron star matter.


    Returns (Implicitly):
    ---------------------
        instance of the class
    """

    EOS_name = None
    tov_mass = None
    EOS_mass_to_rad = None
    EOS_mass_to_bary_mass = None
    threshold_mass = None
    grey_opacity_interp = None
    opacity_data = None

    def __init__(
        self,
        m1=None,
        m2=None,
        c1=None,
        c2=None,
        theta_obs=None,
        Y_e=None,
        m_ej_dyn=None,
        v_ej=None,
        kappa=None,
        m_ej_sec=None,
        m_ej_total=None,
        disk_eff=None,
        EOS=None,
        EOS_path=None,
        kappa_grid_path=None,
        gp_hyperparameter_file=None,
        t=60000.0,
        ra=0.0,
        dec=0.0,
        z=0.001,
        id=None,
        cosmo=cosmos,
        spin1z=None,
        spin2z=None,
        transient_duration=25.0,
        consistency_check=True,
        min_wave=500.0,
        max_wave=12000.0,
        dz_enhancement=1.0,
        thermalisation_eff=0.25,
        mapping_type="coughlin",
        sim_gw=True,
        save=True,
        **kwargs,
    ):
        """Init SAEE viewing-angle class."""
        if id is None:
            self.id = np.random.randint(0, high=2**31)
        else:
            self.id = int(float(id))
        self.t0 = float(t)
        t_1 = Time(t, format="mjd")
        t_1.format = "gps"
        self.t0_gps = t_1.value
        self.ra = float(ra)
        self.dec = float(dec)
        self.num_params = 12
        self.min_wave = float(min_wave)
        self.max_wave = float(max_wave)
        self.mapping_type = str(mapping_type)
        self.dz_enhancement = float(dz_enhancement)
        self.thermalisation_eff = float(thermalisation_eff)
        self.consistency_check = consistency_check
        self.subtype = "Semi-analytic eigenmode expansion with viewing angle."
        self.save = save

        # Handle setup of EOS dependent mapping objects
        if not self.__class__.EOS_name and EOS:
            load_EOS = True
        elif self.__class__.EOS_name and not EOS:
            warnings.warn(
                "Be aware that you are using a pre-loaded EOS.", category=UserWarning
            )
            load_EOS = False
        elif self.__class__.EOS_name != EOS:
            warnings.warn(
                "Be aware that you are changing to a different EOS.",
                category=UserWarning,
            )
            load_EOS = True
        elif self.__class__.EOS_name == EOS:
            warnings.warn("You have already loaded this EOS.", category=UserWarning)
            load_EOS = False
        else:
            raise Exception("You must specify an EOS.")

        if load_EOS is True:
            if EOS_path is None:
                raise Exception(
                    "You must specify the path to an EOS mass-radius table to load."
                )
            self.__class__.EOS_name = EOS
            E1 = eos.get_EOS_table(EOS=EOS, EOS_path=EOS_path)
            self.__class__.tov_mass = eos.get_max_EOS_mass(E1)
            f1 = eos.get_radius_from_EOS(E1)
            self.__class__.EOS_mass_to_rad = f1
            f2 = eos.get_bary_mass_from_EOS(E1)
            self.__class__.EOS_mass_to_bary_mass = f2
            self.__class__.threshold_mass = eos.calculate_threshold_mass(
                self.__class__.tov_mass, f1
            )
            if kappa_grid_path is None:
                raise Exception(
                    "You must specify path to opacity data to construct the Gaussian process object."
                )
            num_data, gp = mappings.construct_opacity_gaussian_process(
                kappa_grid_path, gp_hyperparameter_file
            )
            self.__class__.grey_opacity_interp = gp
            self.__class__.opacity_data = num_data

        self.transient_duration = float(transient_duration)
        self.EOS_name = self.__class__.EOS_name

        # Instantiate parameter values as none to be filled in later
        # Set the parameter names
        self.param1_name = "m1"
        self.param2_name = "m2"
        self.param3_name = "c1"
        self.param4_name = "c2"
        self.param5_name = "theta_obs"
        self.param6_name = "Y_e"
        self.param7_name = "m_ej_dynamic"
        self.param8_name = "v_ej"
        self.param9_name = "kappa"
        self.param10_name = "m_ej_sec"
        self.param11_name = "m_ej_total"
        self.param12_name = "disk_eff"

        for i in range(self.num_params):
            setattr(self, "param{}".format(i + 1), None)

        if m1 is not None:
            self.param1 = float(m1)
        if m2 is not None:
            self.param2 = float(m2)
        if c1 is not None:
            self.param3 = float(c1)
        if c2 is not None:
            self.param4 = float(c2)
        if theta_obs is not None:
            self.param5 = float(theta_obs)
        if Y_e is not None:
            self.param6 = float(Y_e)
        if m_ej_dyn is not None:
            self.param7 = float(m_ej_dyn)
        if v_ej is not None:
            self.param8 = float(v_ej)
        if kappa is not None:
            self.param9 = float(kappa)
        if m_ej_sec is not None:
            self.param10 = float(m_ej_sec)
        if m_ej_total is not None:
            self.param11 = float(m_ej_total)
        if disk_eff is not None:
            self.param12 = float(disk_eff)

        self.draw_parameters()
        self.make_sed()
        super().__init__(float(z), cosmo, sim_gw)

    def draw_parameters(self):
        """Draw parameters not populated by the user.

        Wrapper function to sample the generating parameter distributions for
        BNS kilonova with viewing angle.

        Returns (Implicitly):
        ---------------------
            self.param1: float
                The gravitational mass of the first neutron star.
            self.param2: float
                The gravitational mass of the second neutron star.
            self.param3: float
                The stellar compactness of the first neutron star.
            self.param4: float
                The stellar compactness of the second neutron star.
            self.param5: float
                The viewing angle of the observer with respect to the merger.
            self.param6: float
                The electron fraction along the line of sight of the observer.
            self.param7: float
                The total ejecta mass of the expanding kilonova material.
            self.param8: float
                The mean ejecta velocity of the expanding kilonova material.
            self.param9: float
                The grey opacity of the material along the line of sight of the
                observer.

        """
        if (None in (self.param1, self.param2)) and (
            None in (self.param7, self.param8)
        ):
            (
                self.param1,
                self.param2,
            ) = population.draw_masses_from_EOS_bounds_with_mass_ratio_cut(
                self.__class__.tov_mass
            )

        if None in (self.param3, self.param4):
            self.param3 = eos.compute_compactnesses_from_EOS(
                self.param1, self.__class__.EOS_mass_to_rad
            )
            self.param4 = eos.compute_compactnesses_from_EOS(
                self.param2, self.__class__.EOS_mass_to_rad
            )

        if None in list([self.param5]):
            self.param5 = population.draw_viewing_angle(inclinations=self.param5)

        if None in list([self.param6]):
            self.param6 = mappings.compute_ye_at_viewing_angle(
                self.param5, self.EOS_name
            )

        if None in (self.param7, self.param8):
            self.param7, self.param8 = mappings.map_to_dynamical_ejecta(
                self.param1,
                self.param3,
                self.param2,
                self.param4,
                self.__class__.EOS_mass_to_bary_mass,
            )
            self.param10, self.param11, self.param12 = mappings.map_to_secular_ejecta(
                self.param1,
                self.param3,
                self.param2,
                self.param4,
                self.param7,
                self.__class__.tov_mass,
            )

        if None in list([self.param9]):
            self.param9 = mappings.map_kne_to_grey_opacity_via_gaussian_process(
                self.param11,
                self.param8,
                self.param6,
                self.__class__.grey_opacity_interp,
                self.__class__.opacity_data,
                grey_opacity=self.param9,
            )

        if self.consistency_check is True:
            self.check_kne_priors()

    def make_sed(self, KNE_parameters=None):
        """Create the SED for this kNe instance.

        Wrapper function to send the selected, and default parameters to the
        fortran library which computes the Kilonova SED evolution.

        Parameters:
        -----------
            KNE_parameters: list
                List of parameters needed to generate the kilonova SED. List format
                input instead of inputting mej, vej, etc. separately. Default=None.

        Returns (Implicitly):
        ---------------------
            self.phase: float array
                The phases in days where the kilonova SEDs are defined.
            self.wave: float array
                The discrete wavelengths at which fluxes are defined for a
                given phase.
            self.flux: float array
                The flux values for each combination of phase and wavelength.
        """
        if KNE_parameters is None:
            KNE_parameters = []
            KNE_parameters.append(600.0 / 86400.0)  # starting time [days] 10 minutes
            KNE_parameters.append(self.transient_duration)  # ending time
            KNE_parameters.append(self.param11)  # total ejecta mass
            KNE_parameters.append(self.param8)  # median ejecta velocity
            KNE_parameters.append(1.3)  # nuclear heating rate exponent
            KNE_parameters.append(self.thermalisation_eff)  # thermalization factor
            KNE_parameters.append(self.dz_enhancement)  # DZ enhancement
            KNE_parameters.append(self.param9)  # The grey opacity
            KNE_parameters.append(150.0)  # Initial temperature [K]
            KNE_parameters.append(self.param6)  # Electron fraction
            # Not reading heating rates from file so feed fortran dummy
            # variables
            KNE_parameters.append(
                True
            )  # Flag to use numerical fit nuclear heating rates
            KNE_parameters.append(False)  # Read heating rates variable
            KNE_parameters.append("dummy string")  # Heating rates file
        self.phase, self.wave, self.flux = mw(
            KNE_parameters, self.min_wave, self.max_wave
        )

    def check_kne_priors(
        self,
        m_upper=0.1,
        m_lower=0.001,
        v_upper=0.4,
        v_lower=0.05,
        kappa_lower=0.1,
        max_iter=10000,
    ):
        """Check consistency of parameters with model boundaries.

        Function to see if the fit functions produce values of the ejecta mass
        and ejecta velocity that are broadly consistent with reasonable physical
        limits on their values. This means no more than a quarter of the total
        binary mass as part of the ejecta or that the velocity of the ejecta
        cannot be ultra-relativistic. These limits are quite loose, but seem
        reasonable.
        """
        self.param11 = (
            None if self.param11 > m_upper or self.param11 < m_lower else self.param11
        )
        self.param8 = (
            None if self.param8 > v_upper or self.param8 < v_lower else self.param8
        )
        self.param9 = None if self.param9 < kappa_lower else self.param9

        it = 0
        while None in (self.param11, self.param8, self.param9):
            self.param12 = None
            if it > max_iter:
                for i in range(self.num_params):
                    setattr(self, "param{}".format(i + 1), None)
            (
                self.param1,
                self.param2,
            ) = population.draw_masses_from_EOS_bounds_with_mass_ratio_cut(
                self.tov_mass, mass1=self.param1, mass2=self.param2
            )

            self.param3 = eos.compute_compactnesses_from_EOS(
                self.param1, self.EOS_mass_to_rad
            )
            self.param4 = eos.compute_compactnesses_from_EOS(
                self.param2, self.EOS_mass_to_rad
            )

            self.param5 = population.draw_viewing_angle(inclinations=self.param5)

            self.param6 = mappings.compute_ye_at_viewing_angle(
                self.param5, self.EOS_name, Ye=self.param6
            )

            self.param7, self.param8 = mappings.map_to_dynamical_ejecta(
                self.param1,
                self.param3,
                self.param2,
                self.param4,
                self.EOS_mass_to_bary_mass,
                mej_dyn=self.param7,
                v_ej=self.param8,
            )
            self.param10, self.param11, self.param12 = mappings.map_to_secular_ejecta(
                self.param1,
                self.param3,
                self.param2,
                self.param4,
                self.param7,
                self.tov_mass,
                disk_effs=self.param12,
                m_sec=self.param10,
                m_tot=self.param11,
            )

            self.param9 = mappings.map_kne_to_grey_opacity_via_gaussian_process(
                self.param11,
                self.param8,
                self.param6,
                self.grey_opacity_interp,
                self.opacity_data,
                grey_opacity=self.param9,
            )
            self.param11 = None if self.param11 > m_upper else self.param11
            self.param11 = None if self.param11 < m_lower else self.param11
            self.param8 = None if self.param8 > v_upper else self.param8
            self.param8 = None if self.param8 < v_lower else self.param8
            self.param9 = None if self.param9 < kappa_lower else self.param9
            it += 1


###############################################################################

# The following functions are for future iterations and are in progress.

###############################################################################


def compute_ye_band_factors(self, n_phi=101):
    inclination = self.param5
    phi_grid = (
        np.sort(
            np.arccos(
                2.0 * np.linspace(start=0.0, stop=1.0, num=n_phi, endpoint=True) - 1.0
            )
        )
        + inclination
        - np.pi / 2.0
    )
    ye = mappings.compute_ye_at_viewing_angle(phi_grid, EOS=self.__class__.EOS_name)
    factor = []
    for i, phi in enumerate(phi_grid):
        if i == 0:
            phi_min = phi
            phi_max = phi + (phi_grid[1] - phi) / 2.0
        elif i == n_phi - 1:
            phi_min = phi - (phi - phi_grid[i - 1]) / 2.0
            phi_max = phi
        else:
            phi_min = phi - (phi - phi_grid[i - 1]) / 2.0
            phi_max = phi + (phi_grid[i + 1] - phi) / 2.0

        F_raw, err = quadrature(compute_fphi, phi_min, phi_max, (inclination))
        F = F_raw / np.pi
        factor.append(F)
    fac_array = np.asarray(factor)
    return ye, fac_array


def compute_fphi(phi, inclination):
    if inclination == np.pi / 2.0:
        theta_min = -np.pi / 2.0
        theta_max = np.pi / 2.0
    elif inclination < np.pi / 2.0 and phi < np.pi / 2.0 - inclination:
        theta_min = -np.pi
        theta_max = np.pi
    elif inclination > np.pi / 2.0 and phi > 3.0 * np.pi / 2.0 - inclination:
        theta_min = -np.pi
        theta_max = np.pi
    else:
        targ1 = np.arccos(np.cos(phi) / np.sin(inclination))
        targ2 = -np.arccos(np.cos(phi) / np.sin(inclination))
        x1 = np.cos(targ1) * np.cos(inclination)
        x2 = np.cos(targ2) * np.cos(inclination)
        y1 = np.sin(targ1)
        y2 = np.sin(targ2)
        theta1 = np.arctan2(y1, x1)
        theta2 = np.arctan2(y2, x2)
        theta_min = np.min([theta1, theta2])
        theta_max = np.max([theta1, theta2])
    return np.sin(inclination) * np.sin(phi) * np.sin(phi) * (
        np.sin(theta_max) - np.sin(theta_min)
    ) + (theta_max - theta_min) * np.cos(inclination) * np.cos(phi)
