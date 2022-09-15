"""Module with classes and methods for constructing individual kilonovae
 instances."""
import pkg_resources
import numpy as np
import astropy.units as units
from astropy.time import Time
from sncosmo import TimeSeriesSource, Model
from bnspopkne.macronovae_wrapper import create_SAEE_SEDs
from bnspopkne import equation_of_state as eos
from bnspopkne import mappings
from bnspopkne import population
from tqdm import tqdm
from multiprocessing import Pool


class Setzer2022_kilonova(object):
    """
    Top-level class for kilonovae transients based on Rosswog, et. al 2017
    semi-analytic model for kilonovae spectral energy distributions. With
    added mappings and population priors to generate kilonovae consistent with
    a broad binary neutron star population.

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
        median_ejecta_velocity: float
            The mean ejecta velocity of the expanding kilonova material.
        grey_opacity: float
            The grey opacity of the material along the line of sight of the
            observer.
        secular_ejecta_mass: float
            The secular ejecta mass of the expanding kilonova material.
        total_ejecta_mass: float
            The total ejecta mass of the expanding kilonova material.


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
        EOS_path=None,
        opacity_data_path=None,
        emulator_path=None,
        id=None,
        only_draw_parameters=False,
        **kwargs,
    ):
        # give each transient an id to differentiate simulations
        if id is None:
            self.id = np.random.randint(0, high=2**31)
        else:
            self.id = int(float(id))
        self.number_of_parameters = 12
        self.min_wave = float(min_wave)
        self.max_wave = float(max_wave)
        self.transient_duration = float(transient_duration)
        # Set default data directories
        if EOS_path is None:
            EOS_path = pkg_resources.resource_filename('bnspopkne', "data/mr_sfho_full_right.csv")
        self.EOS_path = EOS_path
        if emulator_path is None:
            emulator_path = pkg_resources.resource_filename('bnspopkne', "data/paper_kernel_hyperparameters.npy")
        self.emulator_path = emulator_path
        if opacity_data_path is None:
            opacity_data_path = pkg_resources.resource_filename('bnspopkne', "data/paper_opacity_data.csv")
        self.opacity_data_path = opacity_data_path

        # Handle setup of EOS dependent mapping objects
        if not self.__class__.EOS_mass_to_rad:
            print('Should only see this once. EOS')
            E1 = eos.get_EOS_table(EOS_path=self.EOS_path)
            self.__class__.tov_mass = eos.get_max_EOS_mass(E1)
            self.__class__.EOS_mass_to_rad = eos.get_radius_from_EOS(E1)
            self.__class__.threshold_mass = mappings.compute_equation_7(
                self.__class__.tov_mass, self.__class__.EOS_mass_to_rad
            )
        # Handle setup of Gaussian process emulator
        if not self.__class__.grey_opacity_emulator:
            print('Should only see this once. GP')
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
        for i in range(self.number_of_parameters):
            setattr(self, "param{}".format(i + 1), None)

        # parse user-provided parameter values
        if ((mass1 is not None) and (mass1 <= self.__class__.tov_mass)):
            self.param1 = float(mass1)
        elif ((mass1 is not None) and (mass1 > self.__class__.tov_mass)):
            raise Exception(f"The provided mass is not compatible with the EOS.")
        if ((mass2 is not None) and (mass2 <= self.__class__.tov_mass)):
            self.param2 = float(mass2)
        elif ((mass2 is not None) and (mass2 > self.__class__.tov_mass)):
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
            source = TimeSeriesSource(self.phase, self.wave, self.flux, name='SAEE kilonova', version='1.0')
            self.model = Model(source=source)  # add dust effect here.
            self.phase = None
            self.wave = None
            self.flux = None
            self.redshift = 0.0
            self.dist_mpc = 1.0e-5
            self.t0 = 0.0  # Currently always start the model at zero phase
            self.tmax = self.model.maxtime()

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
        self.map_mass_to_compactness()
        self.map_to_kilonova_ejecta()
        self.enforce_emulator_bounds()
        self.emulate_grey_opacity()

    def draw_from_binary_population(self):
        """
        Draw from the population priors on masses and inclination. Wrapper
        to population functions to populate instance parameters with draws from
        the population distributions.

        Returns (Implicitly):
        ---------------------
            self.param1: float or Array
                The mass of the primary (larger) neutron star. [solar masses]
            self.param2: float or Array
                The mass of the secondary neutron star. [solar masses]
            self.param5: float or Array
                The viewing angle of the observer w.r.t. the kilonova. [rads]
            self.param12: float or Array
                The unbinding efficiency of the remnant disk.
        """
        # Draw masses based on what is provided.
        if ((self.param1 is None) and (self.param2 is None)):
            (self.param1, self.param2
             ) = population.draw_masses_from_EOS_bounds_with_mass_ratio_cut(
                self.__class__.tov_mass)
        elif ((self.param1 is None) and (self.param2 is not None)):
            self.param1 = population.draw_mass_from_EOS_bounds(
                min(self.__class__tov_mass, (3.0/2.0)*self.param2),
                m_low=self.param2)
        elif ((self.param1 is not None) and (self.param2 is None)):
            self.param2 = population.draw_mass_from_EOS_bounds(
                self.param1, m_low=max(self.param1*(2.0/3.0), 1.0))
        # If not already set, draw the viewing angles
        if self.param5 is None:
            self.param5 = population.draw_viewing_angle()
        # Disk unbinding efficiency
        if self.param12 is None:
            self.param12 = population.draw_disk_unbinding_efficiency()

    def map_mass_to_compactness(self):
        """
        Wrapper function to the equation of state mappings.
        Using the inverse function interpolated from the given EOS mass-radius
        table, compute the compactness for a given mass.

        Returns (Implicitly):
        ---------------------
            self.param3: float or Array
                The compactness of the primary (larger) mass neutron star.
            self.param4: float or Array
                The compactness of the secondary (smaller) mass neutron star.
        """
        # only compute if not already set
        if self.param3 is None:
            self.param3 = eos.compute_compactnesses_from_EOS(
                self.param1, self.__class__.EOS_mass_to_rad
            )
        # only compute if not already set
        if self.param4 is None:
            self.param4 = eos.compute_compactnesses_from_EOS(
                self.param2, self.__class__.EOS_mass_to_rad
            )

    def map_to_kilonova_ejecta(self):
        """
        Wrapper function to compute all of the kilonova ejecta, apart from
        emulation of the grey opacity given the ejecta parameters.

        Returns (Implicitly):
        ---------------------
            self.param6: float or Array
                The mass-weighted electron fraction composition of the ejecta
                along the line-of-sight of the observer viewing-angle.
            self.param7: float or Array
                The mass of the dynamical ejecta component of the kilonova.
                [solar masses]
            self.param8: float or Array
                The median velocity of the dynamical ejecta. [c]
            self.param10: float or Array
                The mass of the secular ejecta component of the kilonova.
                [solar masses]
            self.param11: float or Array
                The total ejecta mass contributing to the kilonova.
                [solar masses]
        """
        # Electron fraction
        if self.param6 is None:
            self.param6 = mappings.compute_equation_10(self.param5)
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

        # Secular Ejecta and Total Ejecta Mass
        total_binary_mass = self.param1 + self.param2
        if self.param10 is None and self.param11 is None:
            self.param10 = mappings.compute_secular_ejecta_mass(total_binary_mass, self.__class__.threshold_mass, self.param12)
            self.param11 = mappings.compute_equation_9(self.param7, self.param10)
        elif self.param10 is None and self.param11 is not None:
            self.param10 = self.param11 - self.param7
            self.param12 = self.param10/mappings.compute_equation_6(total_binary_mass, self.__class__.threshold_mass)
        elif self.param10 is not None and self.param11 is None:
            self.param12 = self.param10/mappings.compute_equation_6(total_binary_mass, self.__class__.threshold_mass)
            self.param11 = mappings.compute_equation_9(self.param7, self.param10)

    def emulate_grey_opacity(self):
        """
        Wrapper function to emulate the grey opacity as the last step in
        populating the transient instance with kilonova parameters.

        Returns (Implicitly):
        ---------------------
            self.param9: float or Array
                The grey opacity of the ejecta.
        """
        if self.param9 is None:
            self.param9 = mappings.emulate_grey_opacity_from_kilonova_ejecta(
                self.param11,
                self.param8,
                self.param6,
                self.__class__.grey_opacity_emulator,
                self.__class__.opacity_data,
            )

    def enforce_emulator_bounds(self):
        """
        Imposes the boundary conditions we are limited to by our training data
        for the Gaussian process emulator.

        If the drawn parameters not fall with the bounds re-draw consistently
        from the population.
        """
        # impose velocity region boundary conditions for heating rates and
        # emulator
        while ((self.param8 < 0.05) or (self.param8 > 0.4) or (self.param11 > 0.08) or (self.param11 < 0.002)):
            self.param1 = None
            self.param2 = None
            self.param3 = None
            self.param4 = None
            self.param7 = None
            self.param8 = None
            self.param10 = None
            # as the disk unbinding efficiency is independent only redraw if
            # the parameters it impacts are outside the bounds
            if ((self.param11 > 0.08) or (self.param11 < 0.002)):
                self.param12 = None
            self.param11 = None
            self.draw_from_binary_population()
            self.map_mass_to_compactness()
            self.map_to_kilonova_ejecta()

    def create_spectral_timeseries(self, KNE_parameters=None):
        """Create the spectral energy density timeseries for this kilonova
        instance.

        First layer wrapper function to FORTRAN library that computes the
        kilonova SED evolution.

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
        self.phase, self.wave, self.flux = create_SAEE_SEDs(
            KNE_parameters, self.min_wave, self.max_wave
        )


class Setzer2022_population_parameter_distribution(object):
    """
    Class to construct population from Setzer et al. 2022.

    This class holds attributes which are vectors of each parameter of the
    population.

    The user can generate a new population, based on size,
    otherwise the population from the paper is imported.

    Returns:
    --------
        Setzer2022_population_parameter_distribution: class instance
            Attributes of the class object are self.paramsi where i is in
            range(12) and are arrays of the distributions corresponding to each
            parameter.

    """
    def __init__(
        self,
        population_size=50000,
        only_draw_parameters=True,
        chunksize=500,
    ):
        """
        Parameters:
        -----------
            population_size: int
                The size of the population to be drawn from the population.
            only_draw_parameters: boolean
                Flag whether to generate only the parameter distributions, or
                also the properties derived from the simulated lightcurves.
                Default is True.
        """
        self.population_size = population_size
        self.number_of_parameters = 12
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
        self.peak_time = np.empty(self.population_size)
        self.peak_abs_lssti = np.empty(self.population_size)
        self.one_mag_peak_time_lssti = np.empty(self.population_size)
        for i in range(self.number_of_parameters):
            setattr(self, "param{}".format(i + 1),
                    np.empty(self.population_size))

        if only_draw_parameters:
            for i in tqdm(range(self.population_size)):
                kilonova = Setzer2022_kilonova(only_draw_parameters=only_draw_parameters)
                for k in range(self.number_of_parameters):
                    getattr(self, f"param{k + 1}")[i] = (
                                            getattr(kilonova, f"param{k + 1}"))
        else:
            with Pool() as p:
                with tqdm(total=int(self.population_size)) as progress_bar:
                    for _ in p.imap_unordered(self.compute_per_kilonova,
                                              list(range(self.population_size)),
                                              chunksize=chunksize):
                        progress_bar.update()

    def compute_per_kilonova(self, id, **kwargs):
        """
        Wrapper function to execute the per-kilonova population draws and
        computation of light curve properties in parallel.

        Returns (Implicitly):
        ---------------------
            self.peak_time: float or Array
                The time of peak of the given kilonova. [days]
            self.peak_abs_lssti: float or Array
                The peak absolute magnitude of the simulated kilonova in the
                LSST i-band.
            self.one_mag_peak_time_lssti: float or Array
                The time [days], the simulated kilonova spends within one mag
                of its peak brightness.
        """
        kilonova = Setzer2022_kilonova(only_draw_parameters=False)
        times = np.linspace(0.0, kilonova.model.maxtime(), 5001)
        lightcurve_abs_i = kilonova.model.bandmag('lssti', "ab", time=times)
        peak_abs_lssti = kilonova.model.source.peakmag('lssti', 'ab', sampling=0.01)
        one_mag_indices = np.nonzero(lightcurve_abs_i <= peak_abs_lssti + 1)
        one_mag_total_time = times[one_mag_indices][-1] - times[one_mag_indices][0]
        self.peak_time[id] = kilonova.model.source.peakphase('lssti', sampling=0.01)
        self.one_mag_peak_time_lssti[id] = one_mag_total_time
        self.peak_abs_lssti[id] = peak_abs_lssti
        for k in range(self.number_of_parameters):
            getattr(self, f"param{k + 1}")[id] = (
                                    getattr(kilonova, f"param{k + 1}"))
