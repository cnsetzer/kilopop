"""Public module for drawing population intrinsic parameters."""
import numpy as np
from bnspopkne import equation_of_state as eos
from bnspopkne import mappings


class Setzer2022_population(object):
    def __init__(self,
                 EOS='sfho',
                 EOS_path=None,
                 gp_hyperparameter_file=None,
                 kappa_grid_path=None,
                 num_samples=50000,
                 transient_duration=21.0,
                 ):
        """Init SAEE viewing-angle class."""
        self.num_transients = num_samples
        self.pre_dist_params = True # relic flag for use with Astrotog software
        self.num_params = 12
        self.min_wave = 500
        self.max_wave = 1200
        self.mapping_type = "coughlin"
        self.dz_enhancement = 1.0
        self.thermalisation_eff = 0.25
        self.consistency_check = True
        E1 = eos.get_EOS_table(EOS=EOS, EOS_path=EOS_path)
        self.tov_mass = eos.get_max_EOS_mass(E1)
        f1 = eos.get_radius_from_EOS(E1)
        self.EOS_mass_to_rad = f1
        f2 = eos.get_bary_mass_from_EOS(E1)
        self.EOS_mass_to_bary_mass = f2
        self.threshold_mass = eos.calculate_threshold_mass(
            self.tov_mass, f1
        )
        if kappa_grid_path is None:
            raise Exception(
                "You must specify path to opacity data to construct the Gaussian process object."
            )
        num_data, gp = mappings.construct_opacity_gaussian_process(
            kappa_grid_path, gp_hyperparameter_file
        )
        self.grey_opacity_interp = gp
        self.opacity_data = num_data

        self.transient_duration = transient_duration
        self.EOS_name = EOS

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

        (
            self.param1,
            self.param2,
        ) = draw_masses_from_EOS_bounds_with_mass_ratio_cut(self.tov_mass, out_shape=(self.num_transients, 1))

        self.param3 = eos.compute_compactnesses_from_EOS(
            self.param1, self.EOS_mass_to_rad
        )
        self.param4 = eos.compute_compactnesses_from_EOS(
            self.param2, self.EOS_mass_to_rad
        )

        self.param5 = draw_viewing_angle(inclinations=self.param5, out_shape=(self.num_transients, 1))

        self.param6 = mappings.compute_ye_at_viewing_angle(
                self.param5, self.EOS_name
            )

        self.param7, self.param8 = mappings.map_to_dynamical_ejecta(
            self.param1,
            self.param3,
            self.param2,
            self.param4,
            self.EOS_mass_to_bary_mass,
        )
        self.param10, self.param11, self.param12 = mappings.map_to_secular_ejecta(
            self.param1,
            self.param3,
            self.param2,
            self.param4,
            self.param7,
            self.tov_mass,
        )

        self.param9 = mappings.map_kne_to_grey_opacity_via_gaussian_process(
            self.param11,
            self.param8,
            self.param6,
            self.grey_opacity_interp,
            self.opacity_data,
            grey_opacity=self.param9,
        )

        m_upper=0.1
        m_lower=0.001
        v_upper=0.4
        v_lower=0.05
        kappa_lower=0.1

        ind1 = np.argwhere(self.param11 > m_upper)
        ind2 = np.argwhere(self.param11 < m_lower)
        ind3 = np.argwhere(self.param8 > v_upper)
        ind4 = np.argwhere(self.param8 < v_lower)
        ind5 = np.argwhere(self.param9 < kappa_lower)
        ind6 = np.union1d(ind1, ind5)
        minds = np.union1d(ind6, ind2)
        vinds = np.union1d(ind3, ind4)
        all_inds = np.union1d(minds, vinds)

        while all_inds.shape[0] > 0:
            print(all_inds.shape)
            for i in range(self.num_params):

                print(getattr(self, f"param{i+1}").shape)
                getattr(self, "param{}".format(i + 1))[all_inds] = None
            (
                self.param1,
                self.param2,
            ) = draw_masses_from_EOS_bounds_with_mass_ratio_cut(self.tov_mass, out_shape=(self.num_transients, 1))

            self.param3 = eos.compute_compactnesses_from_EOS(
                self.param1, self.EOS_mass_to_rad
            )
            self.param4 = eos.compute_compactnesses_from_EOS(
                self.param2, self.EOS_mass_to_rad
            )

            self.param5 = draw_viewing_angle(inclinations=self.param5, out_shape=(self.num_transients, 1))

            self.param6 = mappings.compute_ye_at_viewing_angle(
                    self.param5, self.EOS_name
                )

            self.param7, self.param8 = mappings.map_to_dynamical_ejecta(
                self.param1,
                self.param3,
                self.param2,
                self.param4,
                self.EOS_mass_to_bary_mass,
            )
            self.param10, self.param11, self.param12 = mappings.map_to_secular_ejecta(
                self.param1,
                self.param3,
                self.param2,
                self.param4,
                self.param7,
                self.tov_mass,
            )

            self.param9 = mappings.map_kne_to_grey_opacity_via_gaussian_process(
                self.param11,
                self.param8,
                self.param6,
                self.grey_opacity_interp,
                self.opacity_data,
                grey_opacity=self.param9,
            )
            ind1 = np.argwhere(self.param11 > m_upper)
            ind2 = np.argwhere(self.param11 < m_lower)
            ind3 = np.argwhere(self.param8 > v_upper)
            ind4 = np.argwhere(self.param8 < v_lower)
            ind5 = np.argwhere(self.param9 < kappa_lower)
            ind6 = np.union1d(ind1, ind5)
            minds = np.union1d(ind6, ind2)
            vinds = np.union1d(ind3, ind4)
            all_inds = np.union1d(minds, vinds)


def draw_viewing_angle(inclinations=None, out_shape=1):
    """Set the observer viewing-angle.

    Function draw the observer angle at which the kilonovae is observed. This is
    drawn assuming uniform distribution of binary orbital plane alignment in the
    Universe and making the equivalence of the observer angle and the polar angle
    of the kNe, due to assumed axisymmetry about the normal axis to the binary
    merger plane.

    Parameters:
    -----------
        inclinations: ndarray (optional)
            Array of exisiting inclinations from which to replace certain values
            if set to None. Useful for redrawing a population if certain generated
            kNe were rejected due to parameter consistency checks.

    Returns:
    --------
        inclinations: float or ndarray
            The observer viewing angle, i.e., inclination, with respect to the
            binary merger plane.
    """
    if inclinations is None:
        o_shape = out_shape
    else:
        ind = np.argwhere(np.isnan(inclinations))
        o_shape = inclinations[ind].shape

    theta_obs = np.arccos(2 * np.random.random_sample(size=o_shape) - 1)  # in radians
    if inclinations is None:
        inclinations = theta_obs
    else:
        inclinations[ind] = theta_obs
    return inclinations


def draw_mass_from_EOS_bounds(max_mass, m_low=1.0, mass=None, out_shape=1):
    """Sample a uniform distribution for the mass.

    Parameters:
    -----------
        max_mass: float
            The maximum mass from which to draw the distribution.
        m_low: float
            The minimum mass of the distribution, default 1.0 solar mass.
        mass: float or nd.array
            The
    """
    if mass is None:
        ind = np.array([1])
        o_shape = out_shape
    else:
        ind = np.argwhere(np.isnan(mass))
        o_shape = mass[ind].shape

    new_mass = np.random.uniform(low=m_low, high=max_mass, size=o_shape)

    if mass is None:
        mass = new_mass
    else:
        mass[ind] = new_mass
    return mass


def draw_masses_from_EOS_bounds_with_mass_ratio_cut(max_mass, m_low=1.0, mass_ratio_cut=2.0 / 3.0, mass1=None, mass2=None, out_shape=1):
    """
    Draw neutron star component mass in the source frame give constraints.

    Draw the neutron star component masses assuming a uniform prior over
    the range of masses, as used by LIGO's BNS compact binary waveform search,
    but with a maximum mass set by the chosen EOS.

    Parameters:
    -----------
        max_mass: float
            The maximum TOV mass for the given EOS.
        m_low: float
            The lower bound on the allowable mass of the neutron star.
        mass_ratio_cut: float
            The mass ratio constraint to apply to the mass sampling.
        mass1: nd.array-like
            The set of already provided masses, in solar masses, or default None.
        mass2: nd.array-like
            Same as mass1.
        out_shape: float or shape array
            The number and array shaping structure for the draws of masses.

    Returns (Implicitly):
    ---------------------
        m1: float
            The source-frame gravitational mass of the first neutron star [solar masses].
        m2: float
            The source-frame gravitational mass of the second neutron star [solar masses].
    """
    m1 = draw_mass_from_EOS_bounds(max_mass, mass=mass1, out_shape=out_shape)
    m2 = draw_mass_from_EOS_bounds(max_mass, mass=mass2, out_shape=out_shape)

    if np.isscalar(m1) is True:
        if m1 < m2:
            m1_exch = m1
            m1 = m2
            m2 = m1_exch
        while m2 / m1 < mass_ratio_cut:
            m1 = draw_mass_from_EOS_bounds(max_mass, mass=mass1, out_shape=out_shape)
            m2 = draw_mass_from_EOS_bounds(max_mass, mass=mass2, out_shape=out_shape)
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
            m1[ridx_1] = draw_mass_from_EOS_bounds(max_mass, out_shape=int_shape)
            m2[ridx_1] = draw_mass_from_EOS_bounds(max_mass, out_shape=int_shape)
            ind3 = np.argwhere(m2 > m1)
            # mass ratio cut
            mass_q = np.divide(m2, m1)
            ind4 = np.argwhere(mass_q < mass_ratio_cut)
            # combine indicies into single set to resample
            ridx_1 = np.union1d(ind3, ind4)

    return m1, m2
