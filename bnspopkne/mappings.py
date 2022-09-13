"""Module which contains mappings between intrinsic population parameters and kilonova parameters."""

import numpy as np
import george
from george.modeling import Model
from pandas import read_csv
from bnspopkne import equation_of_state as eos


def compute_equation_6(M_ns_tot, M_thr):
    """
    
    """
    a = -31.335
    b = -0.9760
    c = 1.0474
    d = 0.05957
    m_disk = np.power(
        10.0, a * (1.0 + b * np.tanh((c - (M_ns_tot / M_thr)) / d))
    )
    if m_disk < 1.0e-3:
        m_disk = 1.0e-3
    return m_disk


def compute_secular_ejecta_mass(total_binary_mass, threshold_mass, disk_unbinding_efficiency):
    disk_mass = compute_equation_6(total_binary_mass, threshold_mass)
    secular_ejecta_mass = disk_unbinding_efficiency*disk_mass
    return secular_ejecta_mass


def compute_equation_4(m1, m2, c1, c2):
    """
    Equation 4 in Setzer et al. 2022.

    Coughlin M. W., Dietrich T., Margalit B., Metzger B. D., 2019,
    MNRAS: Letters, 489, L91

    Parameters:
    -----------
        m1: float
            The mass of the primary (more massive) neutron star.
        m2: float
            The mass of the secondary (less massive) neutron star.
        c1: float
            The compactness of the primary (more massive) neutron star.
        c2: float
            The compactness of the secondary (less massive) neutron star.
    Returns:
    --------
        dynamical_ejecta_mass: float
            The dynamical ejecta mass from the fit.
    """
    a = -0.0719
    b = 0.2116
    d = -2.42
    n = -2.905
    dynamical_ejecta_mass = np.power(
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
    return dynamical_ejecta_mass


def compute_equation_5(m1, m2, c1, c2):
    """
    Equation 5 in Setzer et al. 2022.

    Coughlin M. W., Dietrich T., Margalit B., Metzger B. D., 2019,
    MNRAS: Letters, 489, L91

    Parameters:
    -----------
        m1: float
            The mass of the primary (more massive) neutron star.
        m2: float
            The mass of the secondary (less massive) neutron star.
        c1: float
            The compactness of the primary (more massive) neutron star.
        c2: float
            The compactness of the secondary (less massive) neutron star.
    Returns:
    --------
        vej: float
            The median ejecta velocity from the fit.
    """
    e = -0.3090
    f = -1.879
    g = 0.657
    vej = ((e * m1 * (f * c1 + 1.0)) / (m2)) + (g / 2.0)
    +(((e * m2 * (f * c2 + 1.0)) / (m1)) + (g / 2.0))
    return vej


def map_to_dynamical_ejecta(
    mass1,
    comp1,
    mass2,
    comp2,
    mej_dyn=None,
    v_ej=None,
):
    """
    Wrapper for fit functions from various references: Coughlin et. al 2018 etc.
    to map m1,m2,c1,c2 to mej,vej.

    Parameters:
    -----------
        m1: float
            The mass of the primary (more massive) neutron star.
        m2: float
            The mass of the secondary (less massive) neutron star.
        c1: float
            The compactness of the primary (more massive) neutron star.
        c2: float
            The compactness of the secondary (less massive) neutron star.
    Returns:
    --------
.       mej_dyn: float
            Dynamical ejecta mass from the fit.
        v_ej: float
            The median ejecta velocity of the kilonova dynamical ejecta.
    """
    if mej_dyn is None:
        mej_dyn = compute_equation_4(mass1, mass2, comp1, comp2)
    if v_ej is None:
        v_ej = compute_equation_5(mass1, mass2, comp1, comp2)
    return mej_dyn, v_ej


def compute_equation_10(viewing_angle):
    """
    Equation 5 in Setzer et al. 2022.

    Fit function that determines the average electron fraction of the
    material directly in the line of sight of the observer. Determined
    from performing a least squares fit to the BNS merger data from
    Radice et. al 2018.

    Parameters:
    -----------
        viewing_angle: float or ndarray
            The orientation of the kNe w.r.t the observer.

    Returns:
    --------
        electron_fraction: float or ndarray
            The mass-weighted average electron fraction of the material at the
            given viewing_angle or viewing-angle of the observer.
    """
    a = 0.22704
    b = 0.16147
    electron_fraction = a * (np.cos(viewing_angle) ** 2) + b
    return electron_fraction


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
        hyper_vec = np.array([9.05275106, -3.34210729, -0.43019937, -2.93326251])
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


def map_kne_ejecta_to_grey_opacity_via_gaussian_process(
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
    if grey_opacity is None:
        grey_opacity = kappa
    else:
        grey_opacity[ind] = kappa
    if grey_opacity.ndim == 2:
        return grey_opacity
    elif np.isscalar(m_ej_pred):
        return float(grey_opacity)
    else:
        return np.expand_dims(grey_opacity, axis=1)


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
