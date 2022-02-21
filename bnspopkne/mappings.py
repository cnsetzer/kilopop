

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
