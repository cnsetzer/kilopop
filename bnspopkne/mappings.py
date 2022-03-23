import numpy as np
import george
from george.modeling import Model


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
            m_disk = np.power(10.0, a * (1.0 + b * np.tanh((c - (M_tot / M_thr)) / d)))
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
            kappa[ind2] = np.random.normal(loc=k_mean[ind2], scale=np.sqrt(k_var[ind2]))
            ind2 = np.argwhere(kappa < 0.1)
            kit += 1
        if ind2.shape[0] > 0 and self.threshold_opacity is True:
            kappa[ind2] = 0.1

    if self.param9 is None:
        self.param9 = kappa
    else:
        kappa = kappa.reshape((kappa.size, 1))
        self.param9[ind] = kappa


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
    f = interp1d(mass, radius)
    fp = interp1d(mass, np.gradient(radius, mass, edge_order=2))
    return f, fp


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
    f = interp1d(mass, bary_mass)
    fp = interp1d(mass, np.gradient(bary_mass, mass, edge_order=2))
    return f, fp


def get_max_EOS_mass(self):
    """
    Wrapper function to find the Max TOV mass of the given EOS.

    Returns:
    --------
        max_mass: float
            The maximum TOV mass.
    """
    max_mass = max(self.__class__.EOS_table[0]["grav_mass"].values)
    return max_mass


def map_kne_to_grey_opacity(self):
    """
    This is a wrapper function for whatever I get from Oleg and Stephan.

    Returns (Implicitly):
    ---------------------
        self.param9: float
            The grey opacity of the dynamical ejecta material.
    """
    self.map_ye_kne_to_kappa_via_GP()


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

    theta_obs = np.arccos(2 * np.random.random_sample(size=o_shape) - 1)  # in radians
    if self.param5 is None:
        self.param5 = theta_obs
    else:
        self.param5[ind] = theta_obs


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
        while m2 / m1 < mass_ratio_cut:
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
            m1[ridx_1] = np.random.uniform(low=m_low, high=m_high, size=int_shape)
            m2[ridx_1] = np.random.uniform(low=m_low, high=m_high, size=int_shape)
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
