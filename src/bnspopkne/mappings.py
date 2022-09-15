"""Module which contains mappings between intrinsic population parameters and
 kilonova parameters."""

import numpy as np
import george
from george.modeling import Model
from pandas import read_csv
import dill as pickle


def compute_equation_4(mass1, mass2, compactness1, compactness2):
    """
    Equation 4 in Setzer et al. 2022.

    Coughlin M. W., Dietrich T., Margalit B., Metzger B. D., 2019,
    MNRAS: Letters, 489, L91

    Parameters:
    -----------
        mass1: float
            The mass of the primary (more massive) neutron star.
        mass2: float
            The mass of the secondary (less massive) neutron star.
        compactness1: float
            The compactness of the primary (more massive) neutron star.
        compactness2: float
            The compactness of the secondary (less massive) neutron star.
    Returns:
    --------
        dynamical_ejecta_mass: float
            The dynamical ejecta mass from the fit.
    """
    # fit coefficients
    a = -0.0719
    b = 0.2116
    d = -2.42
    n = -2.905
    dynamical_ejecta_mass = np.power(
        10.0,
        (
            ((a * (1.0 - 2.0 * compactness1) * mass1) / (compactness1))
            + b * mass2 * np.power((mass1 / mass2), n)
            + (d / 2.0)
        )
        + (
            ((a * (1.0 - 2.0 * compactness2) * mass2) / (compactness2))
            + b * mass1 * np.power((mass2 / mass1), n)
            + (d / 2.0)
        ),
    )
    return dynamical_ejecta_mass


def compute_equation_5(mass1, mass2, compactness1, compactness2):
    """
    Equation 5 in Setzer et al. 2022.

    Coughlin M. W., Dietrich T., Margalit B., Metzger B. D., 2019,
    MNRAS: Letters, 489, L91

    Parameters:
    -----------
        mass1: float
            The mass of the primary (more massive) neutron star.
        mass2: float
            The mass of the secondary (less massive) neutron star.
        compactness1: float
            The compactness of the primary (more massive) neutron star.
        compactness2: float
            The compactness of the secondary (less massive) neutron star.
    Returns:
    --------
        median_ejecta_velocity: float
            The median ejecta velocity from the fit.
    """
    # fit coefficients
    e = -0.3090
    f = -1.879
    g = 0.657
    median_ejecta_velocity = (((e * mass1 * (f * compactness1 + 1.0)) / (mass2)) + g/2.0
                              + ((e * mass2 * (f * compactness2 + 1.0)) / (mass1)) + g/2.0)
    return median_ejecta_velocity


def compute_equation_6(total_binary_mass, prompt_collapse_threshold_mass):
    """
    Equation 6 from Setzer et al. 2022.

    Compute the remnant disk mass from CITATION.

    Parameters:
    -----------
        total_binary_mass: float or array
            The total mass [solar masses] of the neutron star binary.
        prompt_collapse_threshold_mass: float
            The prompt collapse threshold mass for the given EOS.
    Returns:
    --------
        remnant_disk_mass: float or array
            The mass [solar masses] of the remnant disk after merger.
    """
    # fit coefficients
    a = -31.335
    b = -0.9760
    c = 1.0474
    d = 0.05957
    remnant_disk_mass = np.power(
                                10.0, a * (1.0 + b * np.tanh(
                                 (c - (total_binary_mass /
                                  prompt_collapse_threshold_mass)) / d))
    )
    if remnant_disk_mass < 1.0e-3:
        remnant_disk_mass = 1.0e-3
    return remnant_disk_mass


def compute_equation_7(tov_mass, EOS_mass_to_rad):
    """
    Function to calculate the prompt collapse threshold mass, given the
    maximum TOV mass and the radius at 1.6 solar mass for the chosen EOS.
    Based on Bauswein et al. 2013.

    Parameters:
    -----------
        tov_mass: float
            The maximum mass from solving the TOV equations for given EOS.
        EOS_mass_to_rad: function
            Interpolator function from provided EOS mass-radius curve.
    Returns:
    --------
        prompt_collpase_mass_threshold: float
            Threshold mass in solar masses for prompt blackhole collapse.
    """
    # fit coefficients
    a = 2.38
    b = 3.606
    radius_1_6_msol = EOS_mass_to_rad(1.6)
    prompt_collpase_mass_threshold = ((a - b *
                                      (tov_mass / radius_1_6_msol)) * tov_mass)
    return prompt_collpase_mass_threshold


def compute_equation_8(disk_unbinding_efficiency, disk_mass):
    """
    Compute the total secular ejecta mass given the mass of the remnant disk
    and the assumed unbinding efficiency of the secular processes.

    Parameters:
    -----------
        disk_unbinding_efficiency: float or array
            The fraction of remnant disk matter unbound and contributing to the
            kilonova ejecta mass. [unitless fraction]
        disk_mass: float or array
            The mass of the remnant disk after merger. [solar masses]

    Returns:
    --------
        secular_ejecta_mass: float or Array
            The mass of the kilonova from secular ejection processes.
            [solar masses]
    """
    secular_ejecta_mass = disk_unbinding_efficiency*disk_mass
    return secular_ejecta_mass


def compute_equation_9(dynamical_ejecta_mass, secular_ejecta_mass):
    """
    Compute the total ejecta mass by simple addition of the ejecta mass from
    the dynamical and secular processes involved in the merger.

    Parameters:
    -----------
        dynamical_ejecta_mass: float or Array
            The ejected mass from dynamical processes. [solar masses]
        secular_ejecta_mass: float or Array
            The ejecta mass from secular processes. [solar masses]
    Returns:
    --------
        total_ejecta_mass: float or Array
            The total amount of mass ejected that contributes to the kilonova
            emission. [solar masses]
    """
    total_ejecta_mass = dynamical_ejecta_mass + secular_ejecta_mass
    return total_ejecta_mass


def compute_equation_10(viewing_angle):
    """
    Equation 10 in Setzer et al. 2022.

    Fit function that determines the average electron fraction of the
    material directly in the line of sight of the observer. Determined
    from performing a least squares fit to the BNS merger data from
    Radice et. al 2018.

    Parameters:
    -----------
        viewing_angle: float or ndarray
            The orientation of the kilonova w.r.t the observer.

    Returns:
    --------
        electron_fraction: float or ndarray
            The mass-weighted average electron fraction of the material at the
            given viewing_angle or viewing-angle of the observer.
    """
    # fit coefficients
    a = 0.22704
    b = 0.16147
    electron_fraction = a * (np.cos(viewing_angle) ** 2) + b
    return electron_fraction


def map_to_dynamical_ejecta(
    mass1,
    comp1,
    mass2,
    comp2,
    dynamical_ejecta_mass=None,
    median_ejecta_velocity=None,
):
    """
    Wrapper for fit functions from various references: Coughlin et. al 2018
    to map mass1, mass2, compactness1, compactness2 to mej,
    median_ejecta_velocity.

    Parameters:
    -----------
        mass1: float
            The mass of the primary (more massive) neutron star.
        mass2: float
            The mass of the secondary (less massive) neutron star.
        compactness1: float
            The compactness of the primary (more massive) neutron star.
        compactness2: float
            The compactness of the secondary (less massive) neutron star.
    Returns:
    --------
.       dynamical_ejecta_mass: float
            Dynamical ejecta mass from the fit.
        median_ejecta_velocity: float
            The median ejecta velocity of the kilonova dynamical ejecta.
    """
    if dynamical_ejecta_mass is None:
        dynamical_ejecta_mass = compute_equation_4(mass1, mass2, comp1, comp2)
    if median_ejecta_velocity is None:
        median_ejecta_velocity = compute_equation_5(mass1, mass2, comp1, comp2)
    return dynamical_ejecta_mass, median_ejecta_velocity


def compute_secular_ejecta_mass(total_binary_mass, threshold_mass,
                                disk_unbinding_efficiency):
    """
    Wrapper function for computing the total secular ejecta mass.

    Parameters:
    -----------
        total_binary_mass: float or array
            The total mass [solar masses] of the neutron star binary.
        threshold_mass: float or array
            The mass threshold [solar mass] for prompt collpase after merger.
        disk_unbinding_efficiency: float or array
            The percentage of matter unbound from the remnant disk contributing
            to the radiating ejecta.

    Returns:
    --------
        secular_ejecta_mass: float or array
            The total mass coming from secular processes, post-merger, to add
            to the total mass comprising the kilonova ejecta.
    """
    disk_mass = compute_equation_6(total_binary_mass, threshold_mass)
    secular_ejecta_mass = compute_equation_8(disk_unbinding_efficiency,
                                             disk_mass)
    return secular_ejecta_mass


def construct_opacity_gaussian_process(opacity_data_path,
                                       hyperparameter_file_path):
    """
    Wrapper function to build the gaussian process instance for obtaining a
    grey opacity value for arbitray values of total ejecta mass, median
    ejecta velocity, and electron fraction.

    Parameters:
    -----------
        opacity_data_path: string
            File path for the data from which the grid
            interpolation is built.
        hyperparameter_file_path: string
            File path for the parameters of the Gaussian process from traiing.

    Returns:
    --------
        opacity_GP: george.GP instance
            Trained interpolator function to map
            (total_ejecta_mass, median_ejecta_velocity, electron_fraction) to
            grey opacity.
        grey_opacity_training_data: nd.array
            The training data corresponding to the grey opacity values fit
            from our weighted chi square optimisation.
    """
    # Import the opacity data, and the hyperparameters from training the GP.
    # FUTURE TODO: Consider loading directly the GP object from pickle
    hyper_parameters = np.load(hyperparameter_file_path)
    opacity_array_from_SuperNu_fits = read_csv(opacity_data_path, index_col=0)
    grey_opacity_training_data = opacity_array_from_SuperNu_fits["kappa"].values
    grey_opacity_training_data_uncertainty = opacity_array_from_SuperNu_fits["sigma_kappa"].values
    total_ejecta_masses = opacity_array_from_SuperNu_fits["m_ej"].values
    median_ejecta_velocities = opacity_array_from_SuperNu_fits["v_ej"].values
    electron_fractions = opacity_array_from_SuperNu_fits["Y_e"].values

    kilonova_ejecta_array = np.empty(shape=(len(total_ejecta_masses), 3))
    kilonova_ejecta_array[:, 0] = total_ejecta_masses
    kilonova_ejecta_array[:, 1] = median_ejecta_velocities
    kilonova_ejecta_array[:, 2] = electron_fractions

    kernel_choice = np.var(grey_opacity_training_data) * george.kernels.Matern52Kernel(
        metric=[0.01, 0.05, 0.05], ndim=3
    )

    opacity_GP = george.GP(mean=tanaka_mean_fixed(), kernel=kernel_choice)
    opacity_GP.compute(kilonova_ejecta_array, grey_opacity_training_data_uncertainty)
    opacity_GP.set_parameter_vector(hyper_parameters)
    return grey_opacity_training_data, opacity_GP


def emulate_grey_opacity_from_kilonova_ejecta(
    total_ejecta_mass, median_ejecta_velocity, electron_fraction,
    grey_opacity_gp, opacity_data,
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
    kilonova_ejecta_predictions = np.empty(shape=(1, 3))
    kilonova_ejecta_predictions[0, 0] = total_ejecta_mass
    kilonova_ejecta_predictions[0, 1] = median_ejecta_velocity
    kilonova_ejecta_predictions[0, 2] = electron_fraction

    grey_opacity_mean, grey_opacity_variance = grey_opacity_gp.predict(
        opacity_data, kilonova_ejecta_predictions, return_var=True, cache=True
    )
    # enforce grey opacity thresholding
    grey_opacity = -1.0
    kit = 0
    while grey_opacity < 0.1 and kit < 10000:
        grey_opacity = np.random.normal(loc=grey_opacity_mean,
                                        scale=np.sqrt(grey_opacity_variance))
        kit += 1
    # currently threshold on 0.1
    if grey_opacity < 0.1:
        grey_opacity = 0.1
    return grey_opacity


class tanaka_mean_fixed(Model):
    """Construct piece-wise mean function model based on Tanaka et al. 2019.

    Mean model class to be used with the Gaussian process model of the opacity
    surface. This is based on the work of Tanaka et. al 2019.
    """
    def get_value(self, kilonova_ejecta_array):
        """Get value function in george format.

        Parameters:
        -----------
            self: class instance [implicit]
                Reference to class instance.
            kilonova_ejecta_array: array
                Array of kilonova ejecta parameter pairs.

        Returns:
        --------
            mean_function_value: array
                The mean function value at the given kilonova parameters.
        """
        mean_function_value = np.zeros((len(kilonova_ejecta_array[:, 0]),))
        # low electron fraction
        mean_function_value[kilonova_ejecta_array[:, 2] <= 0.2] = 25.0
        # transition region (linear interpolation)
        mean_function_value[(kilonova_ejecta_array[:, 2] > 0.2) &
                            (kilonova_ejecta_array[:, 2] < 0.25)] = (((-21.0) / (0.05)) * kilonova_ejecta_array[
                             (kilonova_ejecta_array[:, 2] > 0.2) & (kilonova_ejecta_array[:, 2] < 0.25), 2
                            ] + 109.0)
        # mid electron fraction
        mean_function_value[(kilonova_ejecta_array[:, 2] >= 0.25) &
                            (kilonova_ejecta_array[:, 2] <= 0.35)] = 4.0
        # second transition (linear interpolation)
        mean_function_value[(kilonova_ejecta_array[:, 2] > 0.35) &
                            (kilonova_ejecta_array[:, 2] < 0.4)] = (((-3.0) / (0.05)) * kilonova_ejecta_array[
                             (kilonova_ejecta_array[:, 2] > 0.35) & (kilonova_ejecta_array[:, 2] < 0.4), 2
                            ] + 25.0)
        # high electron fraction opacities
        mean_function_value[kilonova_ejecta_array[:, 2] >= 0.4] = 1.0
        return mean_function_value
