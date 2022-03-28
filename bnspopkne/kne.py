import warnings
from copy import deepcopy
import numpy as np
from astropy.constants import c as speed_of_light_ms
import astropy.units as units
from astropy.time import Time
from sncosmo import TimeSeriesSource, Model, read_griddata_ascii
from astropy.cosmology import Planck18 as cosmos
from astropy.cosmology import z_at_value
from sfdmap import SFDMap as sfd
from extinction import fitzpatrick99 as F99
from extinction import apply
from .macronovae_wrapper import make_rosswog_seds as mw

warnings.filterwarnings("ignore", message="ERFA function")

# Set global module constants.
speed_of_light_kms = speed_of_light_ms.to("km/s").value  # Convert m/s to km/s
# Initialze sfdmap for dust corrections
sfdmap = sfd()


class kilonova(em_transient, compact_binary_inspiral):
    """
    Base class for kilonova transients, groups relevant class methods and attributes.
    """

    def __init__(self):
        self.type = "kne"
        em_transient.__init__(self)
        compact_binary_inspiral.__init__(self)


class em_transient(object):
    """
    Base class for transient instances.

    This groups class methods that should be common to all transient models.
    These include assign the transient a peculiar velocity, redshifting, etc.

    """

    def __init__(self):
        source = TimeSeriesSource(self.phase, self.wave, self.flux)
        self.model = deepcopy(Model(source=source))
        # Use deepcopy to make sure full class is saved as attribute of new class
        # after setting model empty phase, wave, flux attributes
        self.phase = None
        self.wave = None
        self.flux = None

    def put_in_universe(
        self, t, ra, dec, z, pec_vel=None, dl=None, cosmo=cosmos, r_v=3.1, id=None,
    ):
        """
        Function to take transient instance and 'place' it into the simulated
        Universe. It sets spacetime location, ra, dec, t0, redshift, etc.

        Input parameters:
        -----------------

        """
        if id is None:
            self.id = np.random.randint(0, high=2 ** 31)
        else:
            self.id = int(id)
        self.t0 = t
        t_1 = Time(t, format="mjd")
        t_1.format = "gps"
        self.t0_gps = t_1.value
        self.ra = ra
        self.dec = dec
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
        self.tmax = t + self.model.maxtime()
        self.extinct_model(r_v=3.1)

    def redshift(self):
        """
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
        """ "
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

    EOS_name = []
    max_mass = []
    EOS_table = []
    EOS_table2 = []
    EOS_mass_to_rad = []
    EOS_mass_to_bary_mass = []
    EOS_mass_to_rad_prime = []
    transient_duration = []
    grey_opacity_interp = []
    opacity_data = []

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
        gp_interp=None,
        spin1z=None,
        spin2z=None,
        only_load_EOS=None,
        transient_duration=25.0,
        consistency_check=True,
        min_wave=500.0,
        max_wave=12000.0,
        dz_enhancement=1.0,
        thermalisation_eff=0.25,
        mapping_type="kruger",
        threshold_opacity=False,
    ):
        self.gp_interp = gp_interp
        self.min_wave = min_wave
        self.max_wave = max_wave
        self.mapping_type = mapping_type
        self.threshold_opacity = threshold_opacity
        if not self.__class__.EOS_name and EOS:
            self.__class__.EOS_name.append(EOS)
            E1, E2 = self.get_EOS_table(EOS=EOS, EOS_path=EOS_path)
            self.__class__.EOS_table.append(E1)
            self.__class__.EOS_table2.append(E2)
            self.__class__.max_mass.append(self.get_max_EOS_mass())
            f, fp, f2, fp2 = self.get_radius_from_EOS()
            self.__class__.EOS_mass_to_rad.append(f)
            self.__class__.EOS_mass_to_rad2.append(f2)
            f, fp, f2, fp2 = self.get_bary_mass_from_EOS()
            self.__class__.EOS_mass_to_bary_mass.append(f)
            self.__class__.EOS_mass_to_bary_mass2.append(f2)
            self.__class__.transient_duration.append(transient_duration)
            if self.gp_interp is True:
                self.__class__.grey_opacity_interp.append(
                    self.get_kappa_GP(kappa_grid_path)
                )
            else:
                self.__class__.grey_opacity_interp.append(
                    self.get_kappa_interp(kappa_grid_path)
                )
        elif not self.__class__.EOS_name and not EOS:
            print("You must first specify the EOS.")
            print(self.__class__.EOS_name, EOS)
            exit()
        elif self.__class__.EOS_name[0] and not EOS:
            pass
        elif self.__class__.EOS_name[0] != EOS:
            self.__class__.EOS_name[0] = EOS
            E1, E2 = self.get_EOS_table(EOS=EOS, EOS_path=EOS_path)
            self.__class__.EOS_table[0] = E1
            self.__class__.EOS_table2[0] = E2
            self.__class__.max_mass[0] = self.get_max_EOS_mass()
            f, fp, f2, fp2 = self.get_radius_from_EOS()
            self.__class__.EOS_mass_to_rad[0] = f
            self.__class__.EOS_mass_to_rad2[0] = f2
            f, fp, f2, fp2 = self.get_bary_mass_from_EOS()
            self.__class__.EOS_mass_to_bary_mass[0] = f
            self.__class__.EOS_mass_to_bary_mass2[0] = f2
            self.__class__.transient_duration[0] = transient_duration
            if self.gp_interp is True:
                self.__class__.grey_opacity_interp[0] = self.get_kappa_GP(
                    kappa_grid_path
                )
            else:
                self.__class__.grey_opacity_interp[0] = self.get_kappa_interp(
                    kappa_grid_path
                )

        elif self.__class__.EOS_name[0] == EOS:
            pass
        else:
            print("You must specify an EOS.")
            exit()

        if only_load_EOS:
            return
        self.num_params = 12
        self.has_gw = True
        try:
            self.transient_duration = self.__class__.transient_duration[0]
        except IndexError:
            self.transient_duration = transient_duration

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

        self.dz_enhancement = dz_enhancement
        self.thermalisation_eff = thermalisation_eff

        self.draw_parameters(consistency_check=consistency_check)

        self.pre_dist_params = False
        self.make_sed()
        self.subtype = "rosswog semi-analytic with viewing angle"
        super().__init__()

    def draw_parameters(self, consistency_check=True):
        """
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

        # Determine output shape for parameters based on size
        if self.number_of_samples > 1:
            self.out_shape = (self.number_of_samples, 1)
        else:
            self.out_shape = None

        if (None in (self.param1, self.param2)) and (
            None in (self.param7, self.param8)
        ):
            self.draw_masses_from_EOS_bounds()

        if None in (self.param3, self.param4):
            self.compute_compactnesses_from_EOS()

        if None in list([self.param5]):
            self.draw_viewing_angle()

        if None in list([self.param6]):
            self.compute_ye_at_viewing_angle()

        if None in (self.param7, self.param8):
            self.map_binary_to_kne()
            self.calculate_secular_ejecta()

        if None in list([self.param9]):
            self.map_kne_to_grey_opacity()

        if consistency_check is True:
            self.check_kne_priors()

    def make_sed(self, KNE_parameters=None):
        """
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
        self, m_upper=0.1, m_lower=0.001, v_upper=0.4, v_lower=0.05, kappa_lower=0.1
    ):
        """
        Function to see if the fit functions produce values of the ejecta mass
        and ejecta velocity that are broadly consistent with reasonable physical
        limits on their values. This means no more than a quarter of the total
        binary mass as part of the ejecta or that the velocity of the ejecta
        cannot be ultra-relativistic. These limits are quite loose, but seem
        reasonable.
        """
        ind1 = np.argwhere(self.param11 > m_upper)
        ind2 = np.argwhere(self.param11 < m_lower)
        ind3 = np.argwhere(self.param8 > v_upper)
        ind4 = np.argwhere(self.param8 < v_lower)
        if self.threshold_opacity is False:
            ind5 = np.argwhere(self.param9 < kappa_lower)
            ind6 = np.union1d(ind1, ind5)
            ind1 = ind6
        minds = np.union1d(ind1, ind2)
        vinds = np.union1d(ind3, ind4)
        all_inds = np.union1d(minds, vinds)

        while all_inds.shape[0] > 0:
            for i in range(self.num_params):
                getattr(self, "param{}".format(i + 1))[all_inds] = None
            self.draw_masses_from_EOS_bounds()
            self.compute_compactnesses_from_EOS()
            self.draw_viewing_angle()
            self.compute_ye_at_viewing_angle()
            self.map_binary_to_kne()
            self.calculate_secular_ejecta()
            self.map_kne_to_grey_opacity()
            ind1 = np.argwhere(self.param11 > m_upper)
            ind2 = np.argwhere(self.param11 < m_lower)
            ind3 = np.argwhere(self.param8 > v_upper)
            ind4 = np.argwhere(self.param8 < v_lower)
            if self.threshold_opacity is False:
                ind5 = np.argwhere(self.param9 < kappa_lower)
                ind6 = np.union1d(ind1, ind5)
                ind1 = ind6
            minds = np.union1d(ind1, ind2)
            vinds = np.union1d(ind3, ind4)
            all_inds = np.union1d(minds, vinds)


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
    ye = compute_ye_at_arbitrary_angle(phi_grid)
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

        F_raw, err = scipy.integrate.quadrature(
            compute_fphi, phi_min, phi_max, (inclination)
        )
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
