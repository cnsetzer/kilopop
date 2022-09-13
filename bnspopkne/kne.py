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

    def __init__(self, z=None, cosmo=None, dl=None):
        """Init of the em_transient class wrapper."""
        source = TimeSeriesSource(self.phase, self.wave, self.flux, name='SAEE kilonova', version='1.0')
        self.model = Model(source=source)  # add dust effect here.
        self.phase = None
        self.wave = None
        self.flux = None
        if z is not None or dl is not None:
            self.put_in_universe(z, cosmo, dl)
        else:
            self.z = 0.0
            self.dist_mpc = 1.0e-5
        self.tmax = self.observer_merger_time + self.model.maxtime()

    def put_in_universe(self, z=None, cosmo=cosmos, dl=None):
        """Place transient instance into the simulated Universe.

        Set properties related to assumed cosmology, i.e., redshift, and
        luminosity distance.

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
        self.dist_pc = self.dist_mpc * 1000000.0  # convert Mpc to pc
        self.redshift()

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
    # Class variables, saves time cacheing these.
    tov_mass = None
    EOS_mass_to_rad = None
    threshold_mass = None
    grey_opacity_emulator = None
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
        disk_unbinding_efficiency=None,
        transient_duration=25.0,
        min_wave=500.0,
        max_wave=12000.0,
        merger_time=60000.0,
        z=None,
        EOS_path="data_directory_placeholder",
        opacity_data_path="data_directory_placeholder",
        emulator_path="data_directory_placeholder",
        cosmo=None,
        id=None,
        sim_gw=False,
        only_draw_parameters=False,
        **kwargs,
    ):
        """Init SAEE viewing-angle class."""
        if id is None:
            self.id = np.random.randint(0, high=2**31)
        else:
            self.id = int(float(id))
        self.observer_merger_time = float(merger_time)
        observer_merger_time_mjd = Time(merger_time, format="mjd")
        observer_merger_time_mjd.format = "gps"
        self.observer_merger_time_gps = observer_merger_time_mjd.value
        self.num_params = 12
        self.min_wave = float(min_wave)
        self.max_wave = float(max_wave)
        self.transient_duration = float(transient_duration)
        self.EOS_path = EOS_path
        self.emulator_path = emulator_path
        self.opacity_data_path = opacity_data_path

        # Handle setup of EOS dependent mapping objects
        if not self.__class__.EOS_mass_to_rad:
            E1 = eos.get_EOS_table(EOS_path=self.EOS_path)
            self.__class__.tov_mass = eos.get_max_EOS_mass(E1)
            self.__class__.EOS_mass_to_rad = eos.get_radius_from_EOS(E1)
            self.__class__.threshold_mass = eos.calculate_threshold_mass(
                self.__class__.tov_mass, self.__class__.EOS_mass_to_rad
            )
        if not self.__class__.grey_opacity_emulator:
            (
             self.__class__.opacity_data, self.__class__.grey_opacity_emulator
            ) = mappings.construct_opacity_gaussian_process(
                self.opacity_data_path, self.emulator_path
            )

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
        for i in range(self.num_params):
            setattr(self, "param{}".format(i + 1), None)

        # parse user-provided parameter values
        if mass1 is not None and mass1 <= self.__class__.tov_mass:
            self.param1 = float(mass1)
        elif mass1 is not None and mass1 > self.__class__.tov_mass:
            raise Exception(f"The provided mass is not compatible with the EOS.")
        if mass2 is not None and mass2 <= self.__class__.tov_mass:
            self.param2 = float(mass2)
        elif mass2 is not None and mass2 > self.__class__.tov_mass:
            raise Exception(f"The provided mass is not compatible with the EOS.")
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
        if disk_unbinding_efficiency is not None:
            self.param12 = float(disk_unbinding_efficiency)

        self.draw_parameters()
        if not only_draw_parameters:
            self.create_spectral_timeseries()
            super().__init__(float(z), cosmo, sim_gw)

    def draw_from_binary_population(self):
        if self.param1 is None and self.param2 is None:
            (
                self.param1,
                self.param2,
            ) = population.draw_masses_from_EOS_bounds_with_mass_ratio_cut(
                self.__class__.tov_mass,
            )
        elif self.param1 is None and self.param2 is not None:
            (
                self.param1,
                _,
            ) = population.draw_masses_from_EOS_bounds_with_mass_ratio_cut(
                min(self.__class__tov_mass, (3.0/2.0)*self.param2),
            )
        elif self.param1 is not None and self.param2 is None:
            (
                self.param2,
                _,
            ) = population.draw_masses_from_EOS_bounds_with_mass_ratio_cut(
                self.param1, m_low=self.param1*(2.0/3.0)
            )
        if self.param3 is None:
            self.param3 = eos.compute_compactnesses_from_EOS(
                self.param1, self.__class__.EOS_mass_to_rad
            )
        if self.param4 is None:
            self.param4 = eos.compute_compactnesses_from_EOS(
                self.param2, self.__class__.EOS_mass_to_rad
            )
        if self.param5 is None:
            self.param5 = population.draw_viewing_angle()
        if self.param6 is None:
            self.param6 = mappings.compute_equation_10(self.param5)

    def draw_parameters(self):
        """Draw parameters not populated by the user.

        Wrapper function to sample the generating parameter distributions from
        Setzer et al 2022.

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
        self.draw_from_binary_population()
        self.map_to_kilonova_ejecta()
        self.enforce_emulator_bounds()

        if self.param9 is None:
            self.param9 = mappings.map_kne_ejecta_to_grey_opacity_via_gaussian_process(
                self.param11,
                self.param8,
                self.param6,
                self.__class__.grey_opacity_emulator,
                self.__class__.opacity_data,
            )

    def map_to_kilonova_ejecta(self):
        # Dynamical Ejecta
        if None in (self.param7, self.param8):
            dynamical_ejecta_mass, median_ejecta_velocity = mappings.map_to_dynamical_ejecta(
                mass1=self.param1,
                comp1=self.param3,
                mass2=self.param2,
                comp2=self.param4,
            )
            if self.param7 is None:
                self.param7 = dynamical_ejecta_mass
            if self.param8 is None:
                self.param8 = median_ejecta_velocity
        # Disk unbinding efficiency
        if self.param12 is None:
            self.param12 = population.draw_disk_unbinding_efficiency()
        # Secular Ejecta and Total Ejecta Mass
        total_binary_mass = self.param1 + self.param2
        if self.param10 is None and self.param11 is None:
            self.param10 = mappings.compute_secular_ejecta_mass(total_binary_mass, self.__class__.threshold_mass, self.param12)
            self.param11 = self.param7 + self.param10
        elif self.param10 is None and self.param11 is not None:
            self.param10 = self.param11 - self.param7
            self.param12 = self.param10/mappings.compute_equation_6(total_binary_mass, self.__class__.threshold_mass)
        elif self.param10 is not None and self.param11 is None:
            self.param12 = self.param10/mappings.compute_equation_6(total_binary_mass, self.__class__.threshold_mass)
            self.param11 = self.param7 + self.param10

    def enforce_emulator_bounds(self):
        # impose velocity region boundary conditions for heating rates and
        # emulator
        while ((self.param8 < 0.05) or (self.param8 > 0.2) or (self.param11 > 0.08) or (self.param11 < 0.002)):
            self.param1 = None
            self.param2 = None
            self.param3 = None
            self.param4 = None
            self.param7 = None
            self.param8 = None
            self.param10 = None
            self.param11 = None
            if (self.param11 > 0.08) or (self.param11 < 0.002):
                self.param12 = None
            self.draw_from_population()
            self.map_to_kilonova_ejecta()

    def create_spectral_timeseries(self, KNE_parameters=None):
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
            KNE_parameters.append(0.25)  # thermalization factor
            KNE_parameters.append(1.0)  # DZ enhancement
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
