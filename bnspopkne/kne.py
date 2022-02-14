__all__ = ["transient_distribution"]
import os
import re
import sys
import warnings
from copy import deepcopy
import numpy as np
import george
from george.modeling import Model
from scipy.interpolate import interp1d, LinearNDInterpolator
from pandas import read_csv
from astropy.constants import c as speed_of_light_ms
from astropy.time import Time
from sncosmo import TimeSeriesSource, Model, zdist, read_griddata_ascii
import opsimsummary as oss
from .functions import bandflux
from astropy.cosmology import Planck15 as cosmos
from sfdmap import SFDMap as sfd
from extinction import fitzpatrick99 as F99
from extinction import apply
from LSSTmetrics.efficiencyTable import EfficiencyTable as eft
from .macronovae_wrapper import make_rosswog_seds as mw
from pycbc.filter import matched_filter
from pycbc.waveform import get_td_waveform
from pycbc import psd
from pycbc.noise import noise_from_psd
from pycbc.detector import Detector

warnings.filterwarnings("ignore", message="ERFA function")
# import warnings
#
# warnings.filterwarnings("ignore", message="numpy.dtype size changed")
# warnings.filterwarnings("ignore", category=DeprecationWarning)
# Set global module constants.
speed_of_light_kms = speed_of_light_ms.to("km/s").value  # Convert m/s to km/s
# Initialze sfdmap for dust corrections
sfdmap = sfd()


def bandflux(band_throughput, SED_model=None, phase=None, ref_model=None):
    """This is wrapper function to compute either the reference system bandflux
       or the bandflux of a source.

    Parameters:
    ----------
    band_throughput : dict of two np.arrays
        This is a dictionary of the wavelengths and transmission fractions that
        characterize the throughput for the given bandfilter.
    SED_model : :obj:, optional
        The sncosmo model object that represents the spectral energy
        distribution of the simulated source.
    phase : float, optional
        This the phase of transient lifetime that is being calculated.
    ref_model : :obj:, optional
        The source model used to calculated the reference bandflux.

    Returns:
    -------
    bandflux : float
        The computed bandflux that is measured by the instrument for a source.

    """
    band_wave = band_throughput["wavelengths"]
    band_tp = band_throughput["throughput"]
    # Get 'reference' SED
    if ref_model:
        flux_per_wave = ref_model.flux(time=2.0, wave=band_wave)

    # Get SED flux
    if SED_model is not None and phase is not None:
        # For very low i.e. zero registered flux, sncosmo sometimes returns
        # negative values so use absolute value to work around this issue.
        flux_per_wave = abs(SED_model.flux(phase, band_wave))

    # Now integrate the combination of the SED flux and the bandpass
    response_flux = flux_per_wave * band_tp
    bandflux = simps(response_flux, band_wave)
    return np.asscalar(bandflux)


class kilonova(transient):
    """
    Base class for kilonova transients, group relevant class methods.
    """

    def __init__(self):
        self.type = "kne"
        super().__init__()


class transient(object):
    """
    Base class for transient instances.

    This groups class methods that should be common to all transient models.
    These include assign the transient a peculiar velocity, redshifting, etc.

    """

    def __init__(self):

        # self.extend_sed_waves()
        source = TimeSeriesSource(self.phase, self.wave, self.flux)
        self.model = deepcopy(Model(source=source))
        # Use deepcopy to make sure full class is saved as attribute of new class
        # after setting model empty phase, wave, flux attributes
        self.phase = None
        self.wave = None
        self.flux = None

    def put_in_universe(self, id, t, ra, dec, z, cosmo=cosmos, r_v=3.1, pec_vel=None, peak=True, sim_gw=True):
        """
        Function to take transient instance and 'place' it into the simulated
        Universe. It sets spacetime location, ra, dec, t0, redshift, etc.

        Input parameters:
        -----------------

        """
        self.id = int(id)
        self.t0 = t
        t_1 = Time(t, format="mjd")
        t_1.format = "gps"
        self.t0_gps = t_1.value
        self.ra = ra
        self.dec = dec
        self.z = z
        self.dist_mpc = cosmo.luminosity_distance(z).value
        self.peculiar_velocity(pec_vel)
        self.redshift(cosmo)
        self.tmax = t + self.model.maxtime()

        if "has_gw" in self.__dict__.keys() and sim_gw is True:
            self.simulate_inspiral_merger()

        self.extinct_model(r_v=3.1)

        if peak is True:
            self.save_peaks(cosmo)
        return self

    def redshift(self, cosmo):
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
        #lumdist = self.dist_mpc * 1e6  # in pc
        #amp = np.power(10.0 / lumdist, 2)
        #self.model.set(amplitude=amp)

        # Current working around for issue with amplitude...
        mapp = cosmo.distmod(self.z).value + self.model.source_peakmag(
             "lsstz", "ab", sampling=0.1)
        self.model.set_source_peakmag(m=mapp, band="lsstz", magsys="ab", sampling=0.1)

    def extinct_model(self, r_v=3.1):
        phases = np.linspace(self.model.mintime(), self.model.maxtime(), num=1000)
        waves = np.linspace(self.model.minwave(), self.model.maxwave(), num=1000)
        unreddend_fluxes = self.model.flux(phases, waves)
        reddend_fluxes = np.empty_like(unreddend_fluxes)
        uncorr_ebv = sfdmap.ebv(self.ra, self.dec, frame="icrs", unit="radian")
        for i, phase in enumerate(phases):
            reddend_fluxes[i, :] = apply(F99(waves, r_v*uncorr_ebv, r_v=r_v), unreddend_fluxes[i, :])

        source = TimeSeriesSource(phases, waves, reddend_fluxes)
        self.extincted_model = deepcopy(Model(source=source))

    def peculiar_velocity(self, pec_vel=None):
        """ "
        Draw from Gaussian peculiar velocity distribution with width 300km/s
        Reference: Hui and Greene (2006) Also, use transient id to set seed for
        setting the peculiar velocity, this is for reproducibility.
        """
        state = np.random.get_state()
        if self.id < pow(2, 32):
            np.random.seed(seed=self.id)
        else:
            print(
                "For some reason this transient id is > 2^32, it is: {}".format(self.id)
            )
        if pec_vel is None:
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

    def extend_sed_waves(self):
        min_wave = np.min(self.wave)
        max_wave = np.max(self.wave)
        delta_wave = abs(self.wave[0] - self.wave[1])
        extension_flux = 0.0
        new_max = 12000.0
        new_min = 500.0

        # insert new max and min
        if max_wave < new_max:
            for wave in np.arange(
                start=max_wave + delta_wave, stop=new_max + delta_wave, step=delta_wave
            ):
                self.wave = np.insert(self.wave, self.wave.size, wave, axis=0)
                self.flux = np.insert(
                    self.flux, self.flux.shape[1], extension_flux, axis=1
                )

        if min_wave > new_min:
            for wave in np.arange(
                start=min_wave - delta_wave, stop=new_min - delta_wave, step=-delta_wave
            ):
                self.wave = np.insert(self.wave, 0, wave, axis=0)
                self.flux = np.insert(self.flux, 0, extension_flux, axis=1)

    def save_peaks(self, cosmo):
        self.peak_lsstu = self.extincted_model.source_peakmag("lsstu", "ab", sampling=0.1)
        self.peak_lsstg = self.extincted_model.source_peakmag("lsstg", "ab", sampling=0.1)
        self.peak_lsstr = self.extincted_model.source_peakmag("lsstr", "ab", sampling=0.1)
        self.peak_lssti = self.extincted_model.source_peakmag("lssti", "ab", sampling=0.1)
        self.peak_lsstz = self.extincted_model.source_peakmag("lsstz", "ab", sampling=0.1)
        self.peak_lssty = self.extincted_model.source_peakmag("lssty", "ab", sampling=0.1)
        self.peak_abs_lsstu = self.model.source_peakabsmag("lsstu", "ab", sampling=0.1, cosmo=cosmo)
        self.peak_abs_lsstg = self.model.source_peakabsmag("lsstg", "ab", sampling=0.1, cosmo=cosmo)
        self.peak_abs_lsstr = self.model.source_peakabsmag("lsstr", "ab", sampling=0.1, cosmo=cosmo)
        self.peak_abs_lssti = self.model.source_peakabsmag("lssti", "ab", sampling=0.1, cosmo=cosmo)
        self.peak_abs_lsstz = self.model.source_peakabsmag("lsstz", "ab", sampling=0.1, cosmo=cosmo)
        self.peak_abs_lssty = self.model.source_peakabsmag("lssty", "ab", sampling=0.1, cosmo=cosmo)


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
        KNE_parameters: list
            List of parameters needed to generate the kilonova SED. List format
            input instead of inputting mej, vej, etc. separately. Default=None.
        parameter_dist: boolean
            Flag to specify if this should just be a sampling distribution
            of the parameters that define the signals without creating the
            events. Default=False.
        num_samples: float
            The number of transients for which to sample the parameter
            distributions. Default=1.

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
    EOS_mass_to_rad2 = []
    EOS_mass_to_bary_mass2 = []
    EOS_mass_to_rad_prime2 = []
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
        KNE_parameters=None,
        parameter_dist=False,
        num_samples=1,
        opac_debug=None,
        consistency_check=True,
        min_wave=500.0,
        max_wave=12000.0,
        dz_enhancement=1.0,
        thermalisation_eff=0.25,
        mapping_type="kruger",
        threshold_opacity=False,
    ):
        self.gp_interp = gp_interp
        self.parameter_dist = parameter_dist
        self.opac_debug = opac_debug
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
        self.number_of_samples = num_samples
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

        if parameter_dist is True:
            if num_samples > 1:
                self.pre_dist_params = True
            else:
                print(
                    "To generate a parameter distribution you need to supply\
                        a number of samples greater than one."
                )
                exit()
            self.subtype = "semi-analytic eigenmode expansion with viewing angle"
            self.type = "parameter distribution"
        else:
            self.pre_dist_params = False
            self.make_sed(KNE_parameters)
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

    def compute_ye_at_viewing_angle(self):
        """
        Fit function that determines the average electron fraction of the
        material directly in the line of sight of the observer. Determined
        from performing a least squares fit to the BNS merger data from
        Radice et. al 2018.

        Returns (Implicitly):
        ---------------------
            self.param6: float
                The electron fraction of the material.
        """
        if self.param6 is None:
            theta_obs = self.param5
        else:
            ind = np.argwhere(np.isnan(self.param6))
            theta_obs = self.param5[ind]
        if self.__class__.EOS_name[0] == "bhblp":
            ye = 0.20640 * (np.cos(theta_obs) ** 2) + 0.16974
        elif self.__class__.EOS_name[0] == "ls220":
            ye = 0.24599 * (np.cos(theta_obs) ** 2) + 0.15083
        elif self.__class__.EOS_name[0] == "dd2":
            ye = 0.19935 * (np.cos(theta_obs) ** 2) + 0.16013
        elif self.__class__.EOS_name[0] == "sfho":
            ye = 0.22704 * (np.cos(theta_obs) ** 2) + 0.16147
        else:
            print(
                "There is a problem with the given EOS: {}".format(
                    self.__class__.EOS_name[0]
                )
            )
        if self.param6 is None:
            self.param6 = ye
        else:
            self.param6[ind] = ye

    def draw_masses_from_EOS_bounds(self, m_low=1.0, mass_ratio_cut=2.0 / 3.0):
        """
        Draw the neutron star component masses assuming a uniform prior over
        the range of masses, as used by LIGO's BNS compact binary waveform search,
        but with a maximum mass set by the chosen EOS.

        Parameters:
        -----------
            m_low: float
                The lower bound on the allowable mass of the neutron star.

        Returns (Implicitly):
        ---------------------
            self.param1: float
                The gravitational mass of the first neutron star.
            self.param2: float
                The gravitational mass of the second neutron star.
        """
        m_high = self.__class__.max_mass[0]
        if self.param1 is None:
            ind = np.array([1, 1])
            o_shape = self.out_shape
        else:
            ind = np.argwhere(np.isnan(self.param1))
            o_shape = self.param1[ind].shape

        if len(ind) > 0:
            m1 = np.random.uniform(low=m_low, high=m_high, size=o_shape)
        else:
            m1 = self.param1

        if self.param2 is None:
            ind2 = np.array([1, 1])
            o_shape = self.out_shape
        else:
            ind2 = np.argwhere(np.isnan(self.param2))
            o_shape = self.param2[ind2].shape

        if len(ind2) > 0:
            m2 = np.random.uniform(low=m_low, high=m_high, size=o_shape)
        else:
            m2 = self.param2

        if np.isscalar(m1) is True:
            if m1 < m2:
                m1_exch = m1
                m1 = m2
                m2 = m1_exch
            while (m2/m1 < mass_ratio_cut):
                m1 = np.random.uniform(low=m_low, high=m_high, size=o_shape)
                m2 = np.random.uniform(low=m_low, high=m_high, size=o_shape)
                if m1 < m2:
                    m1_exch = m1
                    m1 = m2
                    m2 = m1_exch
        else:
            ind3 = np.argwhere(m2 > m1)
            # mass ratio cut
            mass_q = np.divide(m2, m1)
            ind4 = np.argwhere(mass_q < mass_ratio_cut)
            ridx_1 = np.union1d(ind3, ind4)

            # To obtain uniform sampling in the m1,m2 plane resample if m2 > m1
            while ridx_1.shape[0] > 0:
                int_shape = np.shape(m2[ridx_1])
                m1[ridx_1] = np.random.uniform(
                    low=m_low, high=m_high, size=int_shape
                )
                m2[ridx_1] = np.random.uniform(
                    low=m_low, high=m_high, size=int_shape
                )
                ind3 = np.argwhere(m2 > m1)
                # mass ratio cut
                mass_q = np.divide(m2, m1)
                ind4 = np.argwhere(mass_q < mass_ratio_cut)
                # combine indicies into single set to resample
                ridx_1 = np.union1d(ind3, ind4)

        if self.param1 is None:
            self.param1 = m1
        else:
            self.param1[ind] = m1

        if self.param2 is None:
            self.param2 = m2
        else:
            self.param2[ind2] = m2

    def compute_compactnesses_from_EOS(self, sol_2=False):
        """
        Using the equation for compactness from neutron star mass and radius
        compute the compactnesses for both of the component stars in the system.

        Returns (Implicitly):
        ---------------------
            self.param3: float
                The stellar compactness of the first neutron star.
            self.param4: float
                The stellar compactness of the second neturon star.
        """
        G = 13.271317987 * np.power(10, 10)  # units km, M_sol^-1, (km/s)^2
        c = 299792.458  # units km/s
        if self.param3 is None:
            ind = np.array([1, 1])
        else:
            ind = np.argwhere(np.isnan(self.param3))
        if len(ind) > 0:
            if self.param3 is None:
                m1 = self.param1
            else:
                m1 = self.param1[ind]
            if not sol_2:
                R1 = self.__class__.EOS_mass_to_rad[0](m1)  # units km
            else:
                try:
                    R1 = self.__class__.EOS_mass_to_rad2[0](m1)  # units km
                except ValueError:
                    R1 = self.__class__.EOS_mass_to_rad[0](m1)  # units km
            c1 = (G * m1) / ((c ** 2) * R1)
            if self.param3 is None:
                self.param3 = c1
            else:
                self.param3[ind] = c1
        if self.param4 is None:
            ind = np.array([1, 1])
        else:
            ind = np.argwhere(np.isnan(self.param4))
        if len(ind) > 0:
            if self.param4 is None:
                m2 = self.param2
            else:
                m2 = self.param2[ind]
            if not sol_2:
                R2 = self.__class__.EOS_mass_to_rad[0](m2)  # units km
            else:
                try:
                    R2 = self.__class__.EOS_mass_to_rad2[0](m2)  # units km
                except ValueError:
                    R2 = self.__class__.EOS_mass_to_rad[0](m2)  # units km
            c2 = (G * m2) / ((c ** 2) * R2)
            if self.param4 is None:
                self.param4 = c2
            else:
                self.param4[ind] = c2

    def draw_viewing_angle(self):
        """
        Function draw the observer angle at which the kilonovae is observed. This is
        drawn assuming uniform distribution of binary orbital planes in the Universe
        and making the equivalence of the observer angle and the polar angle of the
        kNe, due to assumed axisymmetry about the normal axis to the binary merger
        plane.

        Returns (Implicitly):
        ---------------------
            self.param5: float
                The observer viewing angle with respect to the binary merger plane.
        """
        if self.param5 is None:
            o_shape = self.out_shape
        else:
            ind = np.argwhere(np.isnan(self.param5))
            o_shape = self.param5[ind].shape

        theta_obs = np.arccos(
            2 * np.random.random_sample(size=o_shape) - 1
        )  # in radians
        if self.param5 is None:
            self.param5 = theta_obs
        else:
            self.param5[ind] = theta_obs

    def check_kne_priors(self, m_upper=0.1, m_lower=0.001, v_upper=0.4, v_lower=0.05, kappa_lower=0.1):
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

    def map_binary_to_kne(self):
        """
        Wrapper for fit functions from Coughlin et. al 2018 to map m1,m2,c1,c2 to
        mej,vej.

        Returns (Implicitly):
        ---------------------
            self.param7: float
                The total ejecta mass of the kilonova dynamical ejecta.
            self.param8: float
                The mean ejecta velocity of the kilonova dynamical ejecta.
        """
        which = self.mapping_type
        if which == "coughlin":
            # Fit params from Coughlin et. al 2018
            a = -0.0719
            b = 0.2116
            d = -2.42
            n = -2.905
            e = -0.3090
            f = -1.879
            g = 0.657

        elif which == "radice":
            alpha = -0.657
            beta = 4.254
            gamma = -32.61
            delta = 5.205
            n = -0.773
            e = -0.287
            f = -3.0
            g = 0.494

        elif which == "kruger":
            a = -9.3335
            b = 114.17
            c = -337.56
            n = 1.5465
            e = -0.3090
            f = -1.879
            g = 0.657

        if self.param7 is None:
            ind = np.array([1, 1])
            m1 = self.param1
            m2 = self.param2
            c1 = self.param3
            c2 = self.param4
        else:
            ind = np.argwhere(np.isnan(self.param7))
            m1 = self.param1[ind]
            m2 = self.param2[ind]
            c1 = self.param3[ind]
            c2 = self.param4[ind]
        if len(ind) > 0:
            if which == "coughlin":
                mej = np.power(
                    10.0,
                    (
                        ((a * (1.0 - 2.0 * c1) * m1) / (c1))
                        + b * m2 * np.power((m1 / m2), n)
                        + (d / 2.0)
                    )
                    + (
                        ((a * (1.0 - 2.0 * c2) * m2) / (c2))
                        + b * m1 * np.power((m2 / m1), n)
                        + (d / 2.0)
                    ),
                )
            elif which == "radice":
                bary_m1 = self.__class__.EOS_mass_to_bary_mass[0](m1)
                bary_m2 = self.__class__.EOS_mass_to_bary_mass[0](m2)
                mej = (1.0e-3) * (
                    (
                        alhpa * np.power(m2 / m1, 1.0 / 3) * ((1.0 - 2.0 * c1) / c1)
                        + beta * np.power(m2 / m1, n)
                        + gamma * (1 - m1 / bary_m1)
                    )
                    * bary_m1
                    + (
                        alhpa * np.power(m1 / m2, 1.0 / 3) * ((1.0 - 2.0 * c2) / c2)
                        + beta * np.power(m1 / m2, n)
                        + gamma * (1 - m2 / bary_m2)
                    )
                    * bary_m2
                    + delta
                )
            elif which == "kruger":
                mej = (1.0e-3) * (
                    ((a / c1) + b * np.power(m2 / m1, n) + c * c1) * m1
                    + ((a / c2) + b * np.power(m1 / m2, n) + c * c2) * m2
                )
            if self.param7 is None:
                self.param7 = mej
            else:
                self.param7[ind] = mej

        if self.param8 is None:
            ind = np.array([1, 1])
            m1 = self.param1
            m2 = self.param2
            c1 = self.param3
            c2 = self.param4
        else:
            if np.isscalar(self.param8) is True:
                ind = []
            else:
                ind = np.argwhere(np.isnan(self.param8))
                m1 = self.param1[ind]
                m2 = self.param2[ind]
                c1 = self.param3[ind]
                c2 = self.param4[ind]
        if len(ind) > 0:
            vej = ((e * m1 * (f * c1 + 1.0)) / (m2)) + (g / 2.0)
            +(((e * m2 * (f * c2 + 1.0)) / (m1)) + (g / 2.0))
            if self.param8 is None:
                self.param8 = vej
            else:
                self.param8[ind] = vej

    def get_EOS_table(self, EOS_path=None, EOS="sfho"):
        """
        Wrapper function to read the given Equation of State mass vs. radius
        diagram as pre-computed with a TOV-solver for use with this program.
        The format of the mass vs. radius data should be (radius, grav_mass,
        baryonic_mass,...)

        Parameters:
        -----------
            EOS_path: str
                The location of the file containing the mass vs. radius data.
            EOS: str
                The name of the equation of state.
            crust: boolean
                Flag indicating if the provided equation state contains a
                outer, stiff, neutron star crust.

        Returns:
        --------
            EOS_data: pandas.dataframe
                Dataframe containing the relationship between mass and radius
                for the given EOS.
        """
        EOS_data = read_csv(EOS_path + "mr_{}_full_right.csv".format(EOS))
        EOS_data2 = read_csv(EOS_path + "mr_{}_full_left.csv".format(EOS))
        return EOS_data, EOS_data2

    def get_radius_from_EOS(self):
        """
        Wrapper function to create the interpolation function from the provided
        EOS data to evaluate the mass vs. radius relation for arbitray mass.

        Returns:
        --------
            f: scipy.interpolate instance
                Function which interpolates mass vs. radius for given EOS.
        """
        mass = self.__class__.EOS_table[0]["grav_mass"]
        radius = self.__class__.EOS_table[0]["radius"]
        mass2 = self.__class__.EOS_table2[0]["grav_mass"]
        radius2 = self.__class__.EOS_table2[0]["radius"]
        f = interp1d(mass, radius)
        f2 = interp1d(mass2, radius2)
        fp = interp1d(mass, np.gradient(radius, mass, edge_order=2))
        fp2 = interp1d(mass2, np.gradient(radius2, mass2, edge_order=2))
        return f, fp, f2, fp2

    def get_bary_mass_from_EOS(self):
        """
        Wrapper function to create the interpolation function from the provided
        EOS data to evaluate the mass vs. radius relation for arbitray mass.

        Returns:
        --------
            f: scipy.interpolate instance
                Function which interpolates mass vs. radius for given EOS.
        """
        mass = self.__class__.EOS_table[0]["grav_mass"]
        bary_mass = self.__class__.EOS_table[0]["bary_mass"]
        mass2 = self.__class__.EOS_table2[0]["grav_mass"]
        bary_mass2 = self.__class__.EOS_table2[0]["bary_mass"]
        f = interp1d(mass, bary_mass)
        f2 = interp1d(mass2, bary_mass2)
        fp = interp1d(mass, np.gradient(bary_mass, mass, edge_order=2))
        fp2 = interp1d(mass2, np.gradient(bary_mass2, mass2, edge_order=2))
        return f, fp, f2, fp2

    def get_max_EOS_mass(self):
        """
        Wrapper function to find the Max TOV mass of the given EOS.

        Returns:
        --------
            max_mass: float
                The maximum TOV mass.
        """
        max_mass = max(self.__class__.EOS_table[0]["grav_mass"].values)
        max_mass2 = max(self.__class__.EOS_table2[0]["grav_mass"].values)
        return max(max_mass, max_mass2)

    def map_kne_to_grey_opacity(self):
        """
        This is a wrapper function for whatever I get from Oleg and Stephan.

        Returns (Implicitly):
        ---------------------
            self.param9: float
                The grey opacity of the dynamical ejecta material.
        """

        if self.opac_debug is not None:
            self.param9 = np.ones(self.out_shape) * self.opac_debug
        else:
            if self.gp_interp is True:
                self.map_ye_kne_to_kappa_via_GP()
            else:
                self.map_ye_kne_to_kappa()
            # self.substitute_ye_mapping_to_kappa()

    def map_ye_kne_to_kappa(self):
        """
        Wrapper funciton to use an interpolation instance or other function to
        calculate the corresponding grey opacity that was fit from simulation
        data to triplets of ejecta mass, ejecta velocity, and electron fraction.

        Returns (Implicitly):
        ---------------------
            self.param9: ndarray
                The grey opacity for the instance set of binary and kilonova
                parameters to generate the kilonovae signal.
        """

        if self.param9 is None:
            ind = np.array([1, 1])
            interp_points = np.hstack((self.param11, self.param8, self.param6))
        else:
            ind = np.argwhere(np.isnan(self.param9))
            interp_points = np.hstack(
                (self.param11[ind], self.param8[ind], self.param6[ind])
            )
        # Interpolated logspace grey opacity
        kappa = self.__class__.grey_opacity_interp[0](interp_points)
        # Convert from logspace opacity to real space.
        kappa = np.power(10.0, kappa)
        kappa = kappa.reshape((kappa.size, 1))

        if self.param9 is None:
            self.param9 = kappa
        else:
            self.param9[ind] = kappa

    def simulate_inspiral_merger(self, appx="TaylorF2", dt=1.0 / 2048.0, f_low=25.0):
        """
        Function to calculate the inspiral signal from the BNS merger.

        Returns:
        --------
            self.gw_signal_plus:
            self.gw_signal_cross:
        """
        gps_sec_year = 31536000.0

        self.polarization = 0.0
        wave_p, wave_c = get_td_waveform(
            approximant=appx,
            mass1=self.param1,
            mass2=self.param2,
            spin1z=0.0,
            spin2z=0.0,
            inclination=self.param5,
            delta_t=dt,
            distance=self.dist_mpc,
            f_lower=f_low,
        )
        # avoid invalid time for far future times in lalsimulation
        # but preserve the RA DEC per year positioning by modulating the time
        # of coalescence by gps time in years
        gps_time_mod_year = (
            self.t0_gps / gps_sec_year - int(self.t0_gps / gps_sec_year)
        ) * gps_sec_year
        self.gps_time_mod_year = gps_time_mod_year

        # There is ringing at the end of the template rather than a cutoff at
        # zero time, i.e., time of merger.
        ringing_crop = wave_p.start_time + wave_p.duration
        wave_p = wave_p.crop(0, ringing_crop)
        wave_c = wave_c.crop(0, ringing_crop)
        # Copy the unshifted waveform to save it as the template for matched-filtering
        self.mf_template = deepcopy(wave_p)
        # Shift merger time to 'time of explosion from detected sources'
        wave_p.start_time += gps_time_mod_year
        wave_c.start_time += gps_time_mod_year
        self.gw_signal_plus = wave_p
        self.gw_signal_cross = wave_c

    def get_kappa_interp(self, csv_loc):
        """
        Wrapper function to build the interpolation instance for obtaining a
        grey opacity value for arbitray values of total ejecta mass, median
        ejecta velocity, and electron fraction.

        Parameters:
        -----------
            csv_loc: string
                absolute file path for the data from which the grid
                interpolation is built

        Returns:
        --------
            k_interp: scipy.interpolate: instance
                Interpolator function to map (m_ej_tot,v_ej,Ye) to grey opacity.
        """
        k_df = read_csv(csv_loc, index_col=0)
        masses = k_df["m_ej"].values
        velocities = k_df["v_ej"].values
        electron_fractions = k_df["Y_e"].values
        x = np.empty(shape=(len(masses), 3))
        x[:, 0] = masses
        x[:, 1] = velocities
        x[:, 2] = electron_fractions
        # Interpolate in logspace
        k_grey = np.log10(k_df["kappa"].values)
        k_interp = LinearNDInterpolator(x, k_grey)
        return k_interp

    def get_kappa_GP(self, csv_loc):
        """
        Wrapper function to build the gaussian process instance for obtaining a
        grey opacity value for arbitray values of total ejecta mass, median
        ejecta velocity, and electron fraction.

        Parameters:
        -----------
            csv_loc: string
                absolute file path for the data from which the grid
                interpolation is built

        Returns:
        --------
            opacity_GP: george.GP: instance
                Trained interpolator function to map (m_ej_tot,v_ej,Ye) to grey opacity.
        """
        opac_dataframe = read_csv(csv_loc, index_col=0)
        self.__class__.opacity_data.append(opac_dataframe["kappa"].values)
        opacity_std = opac_dataframe["sigma_kappa"].values
        masses = opac_dataframe["m_ej"].values
        velocities = opac_dataframe["v_ej"].values
        electron_fractions = opac_dataframe["Y_e"].values
        x = np.empty(shape=(len(masses), 3))
        x[:, 0] = masses
        x[:, 1] = velocities
        x[:, 2] = electron_fractions
        kernel_choice = np.var(
            self.__class__.opacity_data[0]
        ) * george.kernels.Matern52Kernel(metric=[0.01, 0.05, 0.05], ndim=3)
        opacity_GP = george.GP(mean=tanaka_mean_fixed(), kernel=kernel_choice)
        opacity_GP.compute(x, opacity_std)
        opacity_GP.set_parameter_vector(
            np.array([9.05275106, -3.34210729, -0.43019937, -2.93326251])
        )
        return opacity_GP

    def map_ye_kne_to_kappa_via_GP(self):
        """
        Wrapper funciton to use an interpolation instance or other function to
        calculate the corresponding grey opacity that was fit from simulation
        data to triplets of ejecta mass, ejecta velocity, and electron fraction.

        Returns (Implicitly):
        ---------------------
            self.param9: ndarray
                The grey opacity for the instance set of binary and kilonova
                parameters to generate the kilonovae signal.
        """

        if self.param9 is None:
            m_ej_pred = self.param11
            v_ej_pred = self.param8
            Y_e_pred = self.param6
        else:
            ind = np.argwhere(np.isnan(self.param9))
            m_ej_pred = self.param11[ind]
            v_ej_pred = self.param8[ind]
            Y_e_pred = self.param6[ind]
        if np.isscalar(m_ej_pred):
            x_pred = np.empty(shape=(1, 3))
            x_pred[0, 0] = m_ej_pred
            x_pred[0, 1] = v_ej_pred
            x_pred[0, 2] = Y_e_pred
        else:
            x_pred = np.empty(shape=(len(m_ej_pred), 3))
            x_pred[:, 0] = m_ej_pred.flatten()
            x_pred[:, 1] = v_ej_pred.flatten()
            x_pred[:, 2] = Y_e_pred.flatten()

        k_mean, k_var = self.__class__.grey_opacity_interp[0].predict(
            self.__class__.opacity_data[0], x_pred, return_var=True, cache=True
        )

        if np.isscalar(k_mean):
            kappa = -1.0
            ind2 = np.array([1, 1])
            kit = 0
            while kappa < 0.1 and kit < 10000:
                kappa = np.random.normal(loc=k_mean, scale=np.sqrt(k_var))
                kit += 1
            if kappa < 0.1 and self.threshold_opacity is True:
                kappa = 0.1
        else:
            kappa = -1.0 * np.ones(shape=(len(k_mean),))
            ind2 = np.argwhere(kappa < 0.1)
            kit = 0
            while ind2.shape[0] > 0 and kit < 10000:
                kappa[ind2] = np.random.normal(
                    loc=k_mean[ind2], scale=np.sqrt(k_var[ind2])
                )
                ind2 = np.argwhere(kappa < 0.1)
                kit += 1
            if ind2.shape[0] > 0 and self.threshold_opacity is True:
                kappa[ind2] = 0.1

        if self.param9 is None:
            self.param9 = kappa
        else:
            kappa = kappa.reshape((kappa.size, 1))
            self.param9[ind] = kappa

    def calculate_threshold_mass(self):
        """
        Function to calculate the prompt collapse threshold mass, given the
        maximum TOV mass and the radius at 1.6 solar mass for the chosen EOS.

        Returns:
        --------
            M_thr: float
                Threshold mass in solar masses.
        """
        a = 2.38
        b = 3.606
        M_tov = self.__class__.max_mass[0]
        R_16 = self.__class__.EOS_mass_to_rad[0](1.6)
        M_thr = (a - b * (M_tov / R_16)) * M_tov
        return M_thr

    def calculate_secular_ejecta(self):
        """
        Function to compute the additional amount of mass resulting from secular
        channels of mass ejection. This is largely modeled as a fraction of the
        remnant disk mass with a floor of 10e-4 solar masses. We take results
        of recent numerical simulations to estimate this total amount of secular
        mass ejected to be in the range of 10-40% of the disk mass.

        Returns (Implicitly):
        ---------------------
        self.param10: ndarray
            Secular component of the total ejecta mass.
        self.param11: ndarray
            The total ejecta mass, i.e., dynamic and secular.
        """
        which = self.mapping_type

        if self.param12 is None:
            disk_eff = np.random.uniform(0.1, 0.4, size=self.out_shape)
            self.param12 = disk_eff
        else:
            dind = np.argwhere(np.isnan(self.param12[:, 0]))[:, 0]
            self.param12[dind] = np.random.uniform(0.1, 0.4, size=dind[:, None].shape)
            disk_eff = self.param12[dind]

        if which == "coughlin":
            a = -31.335
            b = -0.9760
            c = 1.0474
            d = 0.05957
        elif which == "kruger":
            m1 = self.param1
            c1 = self.param3
            a = -8.1324
            c = 1.4820
            d = 1.7784
        elif which == "radice":
            alpha = 0.084
            beta = 0.127
            gamma = 567.1
            delta = 405.14

        M_thr = self.calculate_threshold_mass()

        if self.param10 is None:
            ind = np.array([1, 1])
        else:
            if np.isscalar(self.param10) is True:
                ind = np.array([])
                if self.param11 is None:
                    self.param11 = np.add(self.param7, self.param10)
            else:
                ind = np.argwhere(np.isnan(self.param10[:, 0]))[:, 0]
                if which == "kruger":
                    m1 = self.param1[ind]
                    c1 = self.param3[ind]

        if ind.shape[0] > 0:
            if self.param10 is None:
                M_tot = np.add(self.param1, self.param2)
            else:
                M_tot = np.add(self.param1[ind], self.param2[ind])
            if which == "coughlin":
                m_disk = np.power(
                    10.0, a * (1.0 + b * np.tanh((c - (M_tot / M_thr)) / d))
                )
                if not np.isscalar(m_disk):
                    disk_ind = np.argwhere(m_disk < 1.0e-3)[:, 0]
                    m_disk[disk_ind] = 1.0e-3
                else:
                    if m_disk < 1.0e-3:
                        m_disk = 1.0e-3
            elif which == "kruger":
                m_disk_intermediate = a * c1 + c
                if not np.isscalar(m_disk_intermediate):
                    disk_ind = np.argwhere(m_disk_intermediate < 5 * 1.0e-4)[:, 0]
                    m_disk_intermediate[disk_ind] = 5 * 1.0e-4
                    m_disk = m1 * np.power(m_disk_intermediate, d)
                else:
                    if m_disk_intermediate < 5 * 1.0e-4:
                        m_disk_intermediate = 5 * 1.0e-4
                        m_disk = m1 * np.power(m_disk_intermediate, d)
            elif which == "radice":
                lambda_tilde = 0.0  # Future todo.
                raise NotImplementedError(
                    "The full Radice et al. 2018 disk formulation is not available. Please use options 'coughlin' or 'kruger'."
                )
                m_disk = alpha + beta * np.tanh((lambda_tilde - gamma) / delta)
                if not np.isscalar(m_disk):
                    disk_ind = np.argwhere(m_disk < 5 * 1.0e-4)[:, 0]
                    m_disk[disk_ind] = 5 * 1.0e-4
                else:
                    if m_disk < 5 * 1.0e-4:
                        m_disk = 5 * 1.0e-4

            m_sec = np.multiply(disk_eff, m_disk)
            if self.param10 is None:
                self.param10 = m_sec
                self.param11 = np.add(self.param7, self.param10)
            else:
                self.param10[ind] = m_sec
                self.param11[ind] = np.add(self.param7[ind], self.param10[ind])

        else:
            pass


class tanaka_mean_fixed(Model):
    """
    Mean model class to be used with the Gaussian process model of the opacity
    surface. This is based on the work of Tanaka et. al 2019.


    """

    def get_value(self, x):
        amp = np.zeros((len(x[:, 0]),))
        amp[x[:, 2] <= 0.2] = 25.0
        amp[(x[:, 2] > 0.2) & (x[:, 2] < 0.25)] = ((-21.0) / (0.05)) * x[
            (x[:, 2] > 0.2) & (x[:, 2] < 0.25), 2
        ] + 109.0
        amp[(x[:, 2] >= 0.25) & (x[:, 2] <= 0.35)] = 4.0
        amp[(x[:, 2] > 0.35) & (x[:, 2] < 0.4)] = ((-3.0) / (0.05)) * x[
            (x[:, 2] > 0.35) & (x[:, 2] < 0.4), 2
        ] + 25.0
        amp[x[:, 2] >= 0.4] = 1.0
        return amp


def compute_ye_at_arbitrary_angle(self, angle):
    """
    Fit function that determines the average electron fraction of the
    material directly in the line of sight of the observer. Determined
    from performing a least squares fit to the BNS merger data from
    Radice et. al 2018.

    Returns (Implicitly):
    ---------------------
        self.param6: float
            The electron fraction of the material.
    """

    if self.__class__.EOS_name[0] == "bhblp":
        ye = 0.20640 * (np.cos(theta_obs) ** 2) + 0.16974
    elif self.__class__.EOS_name[0] == "ls220":
        ye = 0.24599 * (np.cos(theta_obs) ** 2) + 0.15083
    elif self.__class__.EOS_name[0] == "dd2":
        ye = 0.19935 * (np.cos(theta_obs) ** 2) + 0.16013
    elif self.__class__.EOS_name[0] == "sfho":
        ye = 0.22704 * (np.cos(theta_obs) ** 2) + 0.16147
    else:
        print(
            "There is a problem with the given EOS: {}".format(
                self.__class__.EOS_name[0]
            )
        )
    return ye


def compute_ye_band_factors(self, n_phi=101):
    inclination = self.param5
    phi_grid = np.sort(np.arccos(
    2.0 * np.linspace(start=0.0, stop=1.0, num=n_phi, endpoint=True) - 1.0
    )) + inclination - np.pi/2.0
    ye = compute_ye_at_arbitrary_angle(phi_grid)
    factor = []
    for i, phi in enumerate(phi_grid):
        if i == 0:
            phi_min = phi
            phi_max = phi + (phi_grid[1] - phi)/2.0
        elif i == n_phi-1:
            phi_min = phi - (phi - phi_grid[i - 1])/2.0
            phi_max = phi
        else:
            phi_min = phi - (phi - phi_grid[i - 1])/2.0
            phi_max = phi + (phi_grid[i+1] - phi)/2.0

        F_raw, err = scipy.integrate.quadrature(compute_fphi, phi_min, phi_max, (inclination))
        F = F_raw / np.pi
        factor.append(F)
    fac_array = np.asarray(factor)
    return ye, fac_array


def compute_fphi(phi, inclination):
    if inclination == np.pi/2.0:
        theta_min = -np.pi/2.0
        theta_max = np.pi/2.0
    elif inclination < np.pi/2.0 and phi < np.pi/2.0 - inclination:
        theta_min = -np.pi
        theta_max = np.pi
    elif inclination > np.pi/2.0 and phi > 3.0*np.pi/2.0 - inclination:
        theta_min = -np.pi
        theta_max = np.pi
    else:
        targ1 = np.arccos(np.cos(phi)/np.sin(inclination))
        targ2 = -np.arccos(np.cos(phi)/np.sin(inclination))
        x1 = np.cos(targ1)*np.cos(inclination)
        x2 = np.cos(targ2)*np.cos(inclination)
        y1 = np.sin(targ1)
        y2 = np.sin(targ2)
        theta1 = np.arctan2(y1, x1)
        theta2 = np.arctan2(y2, x2)
        theta_min = np.min([theta1, theta2])
        theta_max = np.max([theta1, theta2])
    return np.sin(inclination)*np.sin(phi)*np.sin(phi)*(np.sin(theta_max) - np.sin(theta_min)) + (theta_max - theta_min)*np.cos(inclination)*np.cos(phi)
