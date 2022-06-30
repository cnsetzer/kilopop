"""Module which contains mappings between intrinsic population parameters and kilonova parameters."""

import numpy as np
import george
from george.modeling import Model
from pandas import read_csv
from bnspopkne import equation_of_state as eos


class tanaka_mean_fixed(Model):
    """Construct GP mean model based on Tanaka et al. 2019.

    Mean model class to be used with the Gaussian process model of the opacity
    surface. This is based on the work of Tanaka et. al 2019.
    """

    def get_value(self, x):
        """Get value function in george format."""
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


def map_to_secular_ejecta(
    mass1,
    comp1,
    mass2,
    comp2,
    mej_dyn,
    M_thr,
    disk_effs=None,
    m_sec=None,
    m_tot=None,
    mapping_type="coughlin",
    out_shape=1,
):
    """Compute the non-dynamical timescale components of the ejecta.

    Function to compute the additional amount of mass resulting from secular
    channels of mass ejection. This is largely modeled as a fraction of the
    remnant disk mass with a floor of 10e-4 solar masses. We take results
    of recent numerical simulations to estimate this total amount of secular
    mass ejected to be in the range of 10-40% of the disk mass.

    Parameters:
    -----------


    Returns:
    --------
    m_sec: scalar or ndarray
        Secular component of the total ejecta mass in solar masses.
    m_tot: scalar or ndarray
        The total ejecta mass, i.e., dynamic and secular in solar masses.
    """

    if disk_effs is None:
        disk_effs = np.random.uniform(0.1, 0.4, size=out_shape)
    else:
        dind = np.isnan(disk_effs)
        disk_effs[dind] = np.random.uniform(0.1, 0.4, size=dind.shape)

    if mapping_type == "coughlin":
        a = -31.335
        b = -0.9760
        c = 1.0474
        d = 0.05957
    elif mapping_type == "kruger":
        m1 = mass1
        c1 = comp1
        a = -8.1324
        c = 1.4820
        d = 1.7784
    elif mapping_type == "radice":
        alpha = 0.084
        beta = 0.127
        gamma = 567.1
        delta = 405.14

    if m_sec is None:
        ind = np.array([1, 1])
    else:
        if np.isscalar(m_sec) is True:
            ind = np.array([])
            if m_tot is None:
                m_tot = np.add(mej_dyn, m_sec)
        else:
            ind = np.isnan(m_sec[:, 0])[:, 0]
            if mapping_type == "kruger":
                m1 = mass1[ind]
                c1 = comp1[ind]

    if ind.shape[0] > 0:
        if m_sec is None:
            M_ns_tot = np.add(mass1, mass2)
        else:
            M_ns_tot = np.add(mass1[ind], mass2[ind])
        if mapping_type == "coughlin":
            m_disk = np.power(
                10.0, a * (1.0 + b * np.tanh((c - (M_ns_tot / M_thr)) / d))
            )
            if not np.isscalar(m_disk):
                disk_ind = np.argwhere(m_disk < 1.0e-3)[:, 0]
                m_disk[disk_ind] = 1.0e-3
            else:
                if m_disk < 1.0e-3:
                    m_disk = 1.0e-3
        elif mapping_type == "kruger":
            m_disk_intermediate = a * c1 + c
            if not np.isscalar(m_disk_intermediate):
                disk_ind = np.argwhere(m_disk_intermediate < 5 * 1.0e-4)[:, 0]
                m_disk_intermediate[disk_ind] = 5 * 1.0e-4
                m_disk = m1 * np.power(m_disk_intermediate, d)
            else:
                if m_disk_intermediate < 5 * 1.0e-4:
                    m_disk_intermediate = 5 * 1.0e-4
                    m_disk = m1 * np.power(m_disk_intermediate, d)
        elif mapping_type == "radice":
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

        new_m_sec = np.multiply(disk_effs, m_disk)
        if m_sec is None:
            m_sec = new_m_sec
            m_tot = np.add(mej_dyn, new_m_sec)
        else:
            m_sec[ind] = new_m_sec
            m_tot[ind] = np.add(mej_dyn, new_m_sec[ind])

    else:
        pass
    return m_sec, m_tot, disk_effs


def construct_opacity_gaussian_process(csv_loc, hyperparam_file):
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
    # Import the opacity data, and the hyperparameters from training the GP.
    # FUTURE TODO: Consider loading directly the GP object from pickle

    try:
        hyper_vec = np.load(hyperparam_file)
    except (FileNotFoundError, TypeError) as error:
        print(error)
        hyper_vec = np.array([8.49442447, 0.22154227, -0.34043063, -2.99573227])
    opac_dataframe = read_csv(csv_loc, index_col=0)
    grey_opac_vals = opac_dataframe["kappa"].values
    opacity_std = opac_dataframe["sigma_kappa"].values
    masses = opac_dataframe["m_ej"].values
    velocities = opac_dataframe["v_ej"].values
    electron_fractions = opac_dataframe["Y_e"].values

    x = np.empty(shape=(len(masses), 3))
    x[:, 0] = masses
    x[:, 1] = velocities
    x[:, 2] = electron_fractions

    kernel_choice = np.var(grey_opac_vals) * george.kernels.Matern52Kernel(
        metric=[0.01, 0.05, 0.05], ndim=3
    )

    opacity_GP = george.GP(mean=tanaka_mean_fixed(), kernel=kernel_choice)
    opacity_GP.compute(x, opacity_std)
    opacity_GP.set_parameter_vector(hyper_vec)
    return grey_opac_vals, opacity_GP


def map_kne_to_grey_opacity_via_gaussian_process(
    m_tot, v_ej, Y_e, grey_opacity_gp, opacity_data, grey_opacity=None
):
    """
    Wrapper funciton to use an interpolation instance or other function to
    calculate the corresponding grey opacity that was fit from simulation
    data to triplets of ejecta mass, ejecta velocity, and electron fraction.

    Returns (Implicitly):
    ---------------------
        grey_opacity: ndarray
            The grey opacity for the instance set of binary and kilonova
            parameters to generate the kilonovae signal.
    """

    if grey_opacity is None:
        m_ej_pred = m_tot
        v_ej_pred = v_ej
        Y_e_pred = Y_e
    else:
        ind = np.isnan(grey_opacity)
        m_ej_pred = m_tot[ind]
        v_ej_pred = v_ej[ind]
        Y_e_pred = Y_e[ind]
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

    k_mean, k_var = grey_opacity_gp.predict(
        opacity_data, x_pred, return_var=True, cache=True
    )

    if np.isscalar(k_mean):
        kappa = -1.0
        ind2 = np.array([1, 1])
        kit = 0
        while kappa < 0.1 and kit < 10000:
            kappa = np.random.normal(loc=k_mean, scale=np.sqrt(k_var))
            kit += 1
        if kappa < 0.1:
            kappa = 0.1
    else:
        kappa = -1.0 * np.ones(shape=(len(k_mean),))
        ind2 = np.argwhere(kappa < 0.1)
        kit = 0
        while ind2.shape[0] > 0 and kit < 10000:
            kappa[ind2] = np.random.normal(loc=k_mean[ind2], scale=np.sqrt(k_var[ind2]))
            ind2 = np.argwhere(kappa < 0.1)
            kit += 1
        if ind2.shape[0] > 0:
            kappa[ind2] = 0.1
    print(kappa.size)
    if grey_opacity is None:
        grey_opacity = kappa
    else:
        print(ind.size)
        grey_opacity[ind] = kappa
    return np.expand_dims(grey_opacity, axis=1)


def map_to_dynamical_ejecta(
    mass1,
    comp1,
    mass2,
    comp2,
    bary_mass_mapping,
    mej_dyn=None,
    v_ej=None,
    mapping_type="coughlin",
):
    """
    Wrapper for fit functions from various references: Coughlin et. al 2018 etc.
    to map m1,m2,c1,c2 to mej,vej.

    Parameters:
    -----------


    mapping_type: string
        The choice of reference work from which to use the relations mapping
        between these sets of parameters.

    Returns (Implicitly):
    ---------------------
        self.param7: float
            The total ejecta mass of the kilonova dynamical ejecta.
        v_ej: float
            The mean ejecta velocity of the kilonova dynamical ejecta.


    """
    if mapping_type == "coughlin":
        # Fit params from Coughlin et. al 2018
        a = -0.0719
        b = 0.2116
        d = -2.42
        n = -2.905
        e = -0.3090
        f = -1.879
        g = 0.657

    elif mapping_type == "radice":
        alpha = -0.657
        beta = 4.254
        gamma = -32.61
        delta = 5.205
        n = -0.773
        e = -0.287
        f = -3.0
        g = 0.494

    elif mapping_type == "kruger":
        a = -9.3335
        b = 114.17
        c = -337.56
        n = 1.5465
        e = -0.3090
        f = -1.879
        g = 0.657

    if mej_dyn is None:
        ind = np.array([1, 1])
        m1 = mass1
        m2 = mass2
        c1 = comp1
        c2 = comp2
    else:
        ind = np.isnan(mej_dyn)
        m1 = mass1[ind]
        m2 = mass2[ind]
        c1 = comp1[ind]
        c2 = comp2[ind]
    if len(ind) > 0:
        if mapping_type == "coughlin":
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
        elif mapping_type == "radice":
            bary_m1 = bary_mass_mapping(m1)
            bary_m2 = bary_mass_mapping(m2)
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
        elif mapping_type == "kruger":
            mej = (1.0e-3) * (
                ((a / c1) + b * np.power(m2 / m1, n) + c * c1) * m1
                + ((a / c2) + b * np.power(m1 / m2, n) + c * c2) * m2
            )
        if mej_dyn is None:
            mej_dyn = mej
        else:
            mej_dyn[ind] = mej

    if v_ej is None:
        ind = np.array([1, 1])
        m1 = mass1
        m2 = mass2
        c1 = comp1
        c2 = comp2
    else:
        if np.isscalar(v_ej) is True:
            ind = []
        else:
            ind = np.isnan(v_ej)
            m1 = mass1[ind]
            m2 = mass2[ind]
            c1 = comp1[ind]
            c2 = comp2[ind]
    if len(ind) > 0:
        vej = ((e * m1 * (f * c1 + 1.0)) / (m2)) + (g / 2.0)
        +(((e * m2 * (f * c2 + 1.0)) / (m1)) + (g / 2.0))
        if v_ej is None:
            v_ej = vej
        else:
            v_ej[ind] = vej

    return mej_dyn, v_ej


def compute_ye_at_viewing_angle(inclination, EOS, Ye=None):
    """
    Fit function that determines the average electron fraction of the
    material directly in the line of sight of the observer. Determined
    from performing a least squares fit to the BNS merger data from
    Radice et. al 2018.

    Parameters:
    -----------
        inclination: float or ndarray
            The orientation of the kNe w.r.t the observer.

        EOS: string
            The short-hand name of the equation of state to be used when
            performing the mapping.

        Ye: ndarray (optional)
            Exisiting draws for a population of the mass-weighted electron
            fraction. Useful for redrawing values of a population if certain
            kNe were rejected based on consistency checks.

    Returns:
    --------
        Ye: float or ndarray
            The mass-weighted average electron fraction of the material at the
            given inclination or viewing-angle of the observer.
    """
    if Ye is None:
        theta_obs = inclination
    else:
        ind = np.isnan(Ye)
        theta_obs = inclination[ind]

    if EOS == "bhblp":
        new_ye = 0.20640 * (np.cos(theta_obs) ** 2) + 0.16974
    elif EOS == "ls220":
        new_ye = 0.24599 * (np.cos(theta_obs) ** 2) + 0.15083
    elif EOS == "dd2":
        new_ye = 0.19935 * (np.cos(theta_obs) ** 2) + 0.16013
    elif EOS == "sfho":
        new_ye = 0.22704 * (np.cos(theta_obs) ** 2) + 0.16147
    else:
        print("The specified EOS: {}, is not implemented.".format(EOS))
        raise NotImplementedError

    if Ye is None:
        Ye = new_ye
    else:
        Ye[ind] = new_ye
    return Ye
