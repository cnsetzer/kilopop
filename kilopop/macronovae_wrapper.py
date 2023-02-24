"""Public module which computes the thermal spectra from a given set of
kilonova parameters. Wrapper to the FORTRAN library for the SAEE model. """
import numpy as np
from scipy.interpolate import interp1d
from macronova2py import macronova2py as m2p
from astropy.constants import h, c, sigma_sb, k_B, pc


def create_saee_seds(KNE_parameters,
                     min_wave=500.0,
                     max_wave=12000.0,
                     phases=None,
                     wavelengths=None
                     ):
    """
    Function which prepares the inputs from the python code to the fortran
    library to compute the bolometric luminosity evolution for a kilonova given
    the input generating parameters.

    Parameters
    -----------
    KNE_parameters: list
        List in format to pass to fortran code containing all the parmeter
        inputs needed for the fortran library to generate at kilonova
        time-series evolution.
    min_wave: float (optional)
        Minimum wavelength in Angstroms over which to simulate spectra.
    max_wave: float (optional)
        Maximum wavelength in Angstroms over which to simulate spectra.
    phases: np.array-like (optional)
        phases at which to calculate the spectral evolution.
    wavelengths: np.array-like (optional)
        Wavelengths at which to compute the thermal spectrum.

    Returns
    --------
    create_sed_timeseries: function call
        A function call returning a phase, wavelength, and flux array.
    """
    Nt = 5000  # internal grid of times over which to generate the luminosity
    n = len(KNE_parameters) - 4  # number of base kN parameters.
    MNE_parameters = KNE_parameters[0:n]  # sub select parameters
    func_hrate = KNE_parameters[n]  # flag to use functional heating rate
    func_therm = KNE_parameters[n + 1]  # flag to use functional thermalisation
    read_hrate = KNE_parameters[n + 2]  # flag to read heating rates.
    heating_rates_file = KNE_parameters[n + 3]  # file to read from.
    luminosity_array = m2p.calculate_luminosity(
        n, MNE_parameters, func_hrate, func_therm, read_hrate, heating_rates_file, Nt
    )
    return create_sed_timeseries(luminosity_array,
                          min_wave, max_wave, phases, wavelengths)


def create_sed_timeseries(
    luminosity_array, min_wave=500.0, max_wave=12000.0, phases=None,
    wavelengths=None
):
    """
    Function to compute a planck distribution normalized to sigma*T**4/pi

    Parameters
    -----------
    luminosity_array: nd.array
        Time-series evolution of the bolometric luminosity of the simulated
        kilonova source. For use with scaling the Planck distribution.
        Expected shape (n_time, 4).
    min_wave: float (optional)
        Minimum wavelength in Angstroms over which to simulate spectra.
    max_wave: float (optional)
        Maximum wavelength in Angstroms over which to simulate spectra.
    phases: np.array-like (optional)
        phases at which to calculate the spectral evolution.
    wavelengths: np.array-like (optional)
        Wavelengths at which to compute the thermal spectrum.

    Returns
    --------
    phase_days: array
        Array containing all the phases along the lightcurve evolution at
        which an SED is defined. This is given in days. Shape is (n_time,).
    wavelengths: array
        Values in Angstroms of the covered energy spectrum.
        Shape is (n_wave,).
    flux: array
        Energy per area values for every combination of phase and
        wavelengths characterizing the lightcurve of the source.
        [Ergs/s/cm^2/Angstrom]. Shape is (n_time, n_wave).
    """
    day_in_s = 8.64e4  # one day [s]
    Ang_to_cm = 1.0e-8  # angstrom [cm]
    # extract values from luminosity_array
    phase = luminosity_array[:, 0]  # times [seconds]
    luminosity = luminosity_array[:, 1]  # luminosities [erg/s]
    temperature = luminosity_array[:, 2]  # Temperatures [K]

    # if phases specificed, interpolate solution to those times
    if phases is not None:
        luminosity_interpolator = interp1d(phase, y=luminosity)
        luminosity = luminosity_interpolator(phases * day_in_s)
        temperature_interpolator = interp1d(phase, y=temperature)
        temperature = temperature_interpolator(phases * day_in_s)
        phase_days = phases
    else:
        phase_days = phase / day_in_s
    # if wavelengths not provided, setup wavelength array
    if wavelengths is None:
        wavelengths = np.arange(min_wave, max_wave, 10.0)  # Angstrom
    # distance scaling to 10pc at which spectra is computed
    distance_scaling = 10.0 * pc.cgs.value
    # scaling coefficient for the spectra from luminosity solution
    flux_coefficient = np.divide(luminosity,
                                 (4.0 * (distance_scaling ** 2) * sigma_sb.cgs.value *
                                  np.power(temperature, 4)))
    # output flux f [erg s^-1 cm^-2 Ang.^-1]
    lam_cm = np.multiply(wavelengths, Ang_to_cm)
    # output is ergs / s /cm^3
    planck_values = compute_planck_function(lam_cm, temperature)
    # rescale spectrum with luminosity-derived coefficient
    flux = (planck_values.T * flux_coefficient).T
    # convert to ergs/s /cm^2 /Angstrom
    flux *= Ang_to_cm
    return phase_days, wavelengths, flux


def compute_planck_function(lambda_cm, temperature):
    """
    Planck function of wavelength.

    Parameters
    -----------
    lambda_cm: nd.array
        Wavelengths [cm]. Expected shape is (n_wave,).
    temperature: nd.array
        Time-series of temperatures [K] for spectra generation.
        Expected shape is (n_time,).

    Returns
    --------
    planck_values: nd.array
        The value of the Planck function for a given temp. and wavelength.
        Output shape is (n_time, n_wave).
    """
    # cutoff argument to which we set Planck function to zero
    planck_arg_limit = 100.0
    # utilize linear algebra for fast generation, i.e., broadcasting
    lambda_cm = np.expand_dims(lambda_cm, axis=0)
    temperature = np.expand_dims(temperature, axis=1)
    # pre-compute exponent argument for Planck function
    planck_arg = np.divide(h.cgs.value * c.cgs.value,
                           (k_B.cgs.value * temperature * lambda_cm))
    # construct wavelength array for each temperature
    lam_planck = np.tile(lambda_cm, (planck_arg.shape[0], 1))
    # initialize Planck function array
    planck_values = np.zeros(planck_arg.shape)
    # compute Planck function
    planck_values[planck_arg <= planck_arg_limit] = np.divide(
        (2.0 * h.cgs.value * (c.cgs.value ** 2)) /
        np.power(lam_planck[planck_arg <= planck_arg_limit], 5),
        (np.exp(planck_arg[planck_arg <= planck_arg_limit]) - 1.0),
    )
    return planck_values
