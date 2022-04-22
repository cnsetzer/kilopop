import numpy as np


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
