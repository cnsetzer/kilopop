"""Module with classes and methods for constructing individual kilonovae
 instances."""
import warnings
import numpy as np
from astropy.constants import c as speed_of_light_ms
import astropy.units as units
from astropy.time import Time
from sncosmo import TimeSeriesSource, Model
from astropy.cosmology import Planck18 as cosmos
from astropy.cosmology import z_at_value
from bnspopkne.macronovae_wrapper import make_rosswog_seds as mw
from bnspopkne.inspiral import compact_binary_inspiral
from bnspopkne import equation_of_state as eos
from bnspopkne import mappings
from bnspopkne import population
from scipy.integrate import quadrature

# Filter warnings from Numpy
warnings.filterwarnings("ignore", message="ERFA function")

# Set global module constants.
speed_of_light_kms = speed_of_light_ms.to("km/s").value  # Convert m/s to km/s


class em_transient(object):
    """
    Base class for optical transient instances.

    This groups class methods that should be common to all transient models.
    These include assign the transient a peculiar velocity, redshifting, etc.

    """

    def __init__(self, z=None, cosmo=cosmos, dl=None):
        """Init of the em_transient class wrapper."""
        source = TimeSeriesSource(self.phase, self.wave, self.flux, name='SAEE kilonova', version='1.0')
        self.model = Model(source=source)  # add dust effect here.
        self.phase = None
        self.wave = None
        self.flux = None
        if z is not None:
            self.put_in_universe(z, cosmo, dl)

    def put_in_universe(self, z=None, cosmo=cosmos, dl=None):
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
        self.redshift()
        self.tmax = self.t0 + self.model.maxtime()

    def redshift(self):
        """Redshift the source.

        Wrapper function to redshift the spectrum of the transient instance,
        and scale the flux according to the luminosity distance.

        Returns:
        -------
            self: modified self class instance
        """
        # Note that it is necessary to scale the amplitude relative to the 10pc
        # (i.e. 10^2 in the following eqn.) placement of the SED currently
        amp = np.power(10.0 / self.dist_pc, 2)
        self.model.update({'z': self.z, 'amplitude': amp})

        # Current working around for issue with amplitude...
        # mapp = 5.0 * np.log10(self.dist_pc / 10.0) + self.model.source_peakmag(
        #     "lsstz", "ab", sampling=0.05
        # )
        # self.model.set_source_peakmag(m=mapp, band="lsstz", magsys="ab", sampling=0.05)


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
        mass1: float
            The gravitational mass of the first neutron star.
        mass2: float
            The gravitational mass of the second neutron star.
        compactness1: float
            The stellar compactness of the first neutron star.
        compactness2: float
            The stellar compactness of the second neutron star.
        viewing_angle: float
            The viewing angle of the observer with respect to the merger.
        electron_fraction: float
            The electron fraction along the line of sight of the observer.
        dynamical_ejecta_mass: float
            The dynamic ejecta mass of the expanding kilonova material.
        vej: float
            The mean ejecta velocity of the expanding kilonova material.
        grey_opacity: float
            The grey opacity of the material along the line of sight of the
            observer.
        secular_ejecta_mass: float
            The secular ejecta mass of the expanding kilonova material.
        total_ejecta_mass: float
            The total ejecta mass of the expanding kilonova material.
        EOS: str
            The name of the equation of state of the neutron star matter.


    Returns (Implicitly):
    ---------------------
        Kilonova instance of the class
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
        mass1=None,
        mass2=None,
        compactness1=None,
        compactness2=None,
        viewing_angle=None,
        electron_fraction=None,
        dynamical_ejecta_mass=None,
        median_ejecta_velocity=None,
        grey_opacity=None,
        secular_ejecta_mass=None,
        total_ejecta_mass=None,
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
        self.param12_name = "disk_eff"

        for i in range(self.num_params):
            setattr(self, "param{}".format(i + 1), None)

        if mass1 is not None:
            self.param1 = float(mass1)
        if mass2 is not None:
            self.param2 = float(mass2)
        if compactness1 is not None:
            self.param3 = float(compactness1)
        if compactness2 is not None:
            self.param4 = float(compactness2)
        if viewing_angle is not None:
            self.param5 = float(viewing_angle)
        if electron_fraction is not None:
            self.param6 = float(electron_fraction)
        if dynamical_ejecta_mass is not None:
            self.param7 = float(dynamical_ejecta_mass)
        if median_ejecta_velocity is not None:
            self.param8 = float(median_ejecta_velocity)
        if grey_opacity is not None:
            self.param9 = float(grey_opacity)
        if secular_ejecta_mass is not None:
            self.param10 = float(secular_ejecta_mass)
        if total_ejecta_mass is not None:
            self.param11 = float(total_ejecta_mass)
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
            KNE_parameters.append(600.0 / 86400.0)  # starting time [days] (10 minutes)
            KNE_parameters.append(self.transient_duration)  # ending time [days]
            KNE_parameters.append(self.param11)  # total ejecta mass
            KNE_parameters.append(self.param8)  # median ejecta velocity
            KNE_parameters.append(1.3)  # heating rate exponent (not used)
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
            KNE_parameters.append("placeholder string")  # Heating rates file
        self.phase, self.wave, self.flux = mw(
            KNE_parameters, self.min_wave, self.max_wave
        )

    def check_kne_priors(
        self,
        m_upper=0.1,
        m_lower=0.001,
        v_upper=0.4,
        v_lower=0.05,
        opacity_lower_bound=0.1,
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
        self.param9 = None if self.param9 < opacity_lower_bound else self.param9

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
                median_ejecta_velocity=self.param8,
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
            self.param9 = None if self.param9 < opacity_lower_bound else self.param9
            it += 1
