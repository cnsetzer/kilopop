import numpy as np
from scipy.interpolate import interp1d
from macronova2py import macronova2py as m2p

#####################
# physics constants #
#####################
parsec = 3.08567758e18  # parsec [cm]
clight = 2.99792458e10  # speed of light [cm/s]
sigma = 5.670374419e-5  # Stefan-Boltzmann constant [erg/(s*cm2*K4)]
day_in_s = 8.64e4  # one day [s]
Ang_to_cm = 1.0e-8  # angstrom [cm]
Robs = 10.0 * parsec  # distance for absolute magnitudes
hplanck = 6.62607015e-27  # Planck constant [erg*s]
kB = 1.38064852e-16  # Boltzmann constant [erg/K]


def make_rosswog_seds(KNE_parameters,
                      min_wave=500.0,
                      max_wave=12000.0,
                      phases=None,
                      wavelengths=None
                      ):
    """
    Function which prepares the inputs from the python code to the fortran
    library to compute the bolometric luminosity evolution for a kilonova given
    the input generating parameters.

    Parameters:
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
    Returns:
    --------
        sed_timeseries: function call
            A function call returning a phase, wavelength, and flux array.
    """
    Nt = 5000  # internal grid of times over which to generate the luminosity
    n = len(KNE_parameters) - 3  # number of parameters for base kN params.
    MNE_parameters = KNE_parameters[0:n]  # sub select parameters
    func_hrate = KNE_parameters[n]  # flag to use functional heating rate
    read_hrate = KNE_parameters[n + 1]  # flag to read heating rates.
    heating_rates_file = KNE_parameters[n + 2]  # file to read from.
    luminosity = m2p.calculate_luminosity(
        n, MNE_parameters, func_hrate, read_hrate, heating_rates_file, Nt
    )
    return sed_timeseries(luminosity, min_wave, max_wave, phases, wavelengths)


def sed_timeseries(
    luminosity, min_wave=500.0, max_wave=12000.0, phases=None, wavelengths=None
):
    """
    Function to compute a planck distribution normalized to sigma*T**4/pi

    Parameters:
    -----------
        luminosity: nd.array
            Time-series evolution of the bolometric luminosity of the simulated
            kilonova source. For use with scaling the Planck distribution. Expected shape (n_time, 4).
        min_wave: float (optional)
            Minimum wavelength in Angstroms over which to simulate spectra.
        max_wave: float (optional)
            Maximum wavelength in Angstroms over which to simulate spectra.
        phases: np.array-like (optional)
            phases at which to calculate the spectral evolution.
        wavelengths: np.array-like (optional)
            Wavelengths at which to compute the thermal spectrum.

    Returns:
    --------
        ti: array
            Array containing all the phases along the lightcurve evolution at
            which an SED is defined. This is given in seconds. Shape is (n_time,).
        wavelengths: array
            Values in Angstroms of the covered energy spectrum. Shape is (n_wave,).
        flux: array
            Energy per area values for every combination of phase and wavelengths
            characterizing the lightcurve of the source, [Ergs/s/cm^2/Angstrom]. Shape is (n_time, n_wave).
    """
    # extract values from luminosity
    ti = luminosity[:, 0]  # times [seconds]
    Li = luminosity[:, 1]  # luminosities [erg/s]
    Ti = luminosity[:, 2]  # Temperatures [K]

    if phases is not None:
        Linterp = interp1d(ti, y=Li)
        Li = Linterp(phases * day_in_s)
        Tinterp = interp1d(ti, y=Ti)
        Ti = Tinterp(phases * day_in_s)
        ti_days = phases
    else:
        ti_days = ti / day_in_s
    # if wavelengths not provided, setup wavelength array
    if wavelengths is None:
        wavelengths = np.arange(min_wave, max_wave, 10.0)  # Angstrom
    # output array,
    Coef = np.divide(Li, (4.0 * (Robs ** 2) * sigma * np.power(Ti, 4)))
    print(Coef.shape)
    # output flux f [erg s^-1 cm^-2 Ang.^-1]
    lam_cm = np.multiply(wavelengths, Ang_to_cm)
    # use normalisation coefficient and integrate the per steraidan blam over the outward solid angle
    # i.e., that which is toward the observer
    blam_test = blam(lam_cm, Ti)
    print(blam_test.shape)
    flux = (blam_test.T * Coef).T
    flux *= Ang_to_cm
    return ti_days, wavelengths, flux


def blam(lam, T):
    """
    Planck function of wavelength.

    Parameters:
    -----------
        lam: np.array
            Wavelengths in cm. Expected shape is (n_wave,).
        T: np.array
            Time-series of temperatures for model evoluation. Expected shape is (n_time,).

    Returns:
    --------
        Planck: np.array
            The value of the Planck function for a given temp. and wavelength. Output shape is (n_time, n_wave).
    """
    # cutoff argument to which we set Planck function to zero
    x_cut = 100.0
    lam = np.expand_dims(lam, axis=0)
    T = np.expand_dims(T, axis=1)
    x = np.divide(hplanck * clight, (kB * T * lam))
    lam_planck = np.tile(lam, (x.shape[0], 1))
    Planck = np.zeros(x.shape)
    Planck[x <= x_cut] = np.divide(
        (2.0 * hplanck * clight ** 2) / np.power(lam_planck[x <= x_cut], 5),
        (np.exp(x[x <= x_cut]) - 1.0),
    )
    return Planck
