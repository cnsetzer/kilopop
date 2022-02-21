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
