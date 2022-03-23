class compact_binary_inspiral(object):
    def __init__(self, appx="TaylorF2", dt=1.0 / 2048.0, f_low=25.0, **kwargs):
        self.appx = appx
        self.dt = dt
        self.f_low = f_low
        self.simulate_inspiral_merger()
        super.__init__()

    def simulate_inspiral_merger(self):
        """
        Function to calculate the inspiral signal from the BNS merger.

        Parameters:
        -----------
            appx: string
                Name of approximant from pycbc to use as the gw template.
            dt: float
                Sampling rate delta time for waveform generation.
            f_low: float
                Lower frequency at which to being the waveform.

        Generates:
        --------
            self.gw_signal_plus: Pycbc TimeSeries object
                The plus polarization of the (2,2) mode gravitational wave.
            self.gw_signal_cross: Pycbc TimeSeries object
                The cross polarization of the (2,2) mode gravitational wave.
        """
        gps_sec_year = 31536000.0

        self.polarization = 0.0
        wave_p, wave_c = get_td_waveform(
            approximant=self.appx,
            mass1=self.param1,
            mass2=self.param2,
            spin1z=0.0,
            spin2z=0.0,
            inclination=self.param5,
            delta_t=self.dt,
            distance=self.dist_mpc,
            f_lower=self.f_low,
        )
        # avoid invalid time for far future times in lalsimulation
        # but preserve the RA DEC per year positioning by modulating the time
        # of coalescence by gps time in years
        gps_time_mod_year = (
            self.t0_gps / gps_sec_year - int(self.t0_gps / gps_sec_year)
        ) * gps_sec_year
        self.gps_time_mod_year = gps_time_mod_year

        # There is ringing at the end of the template rather than a cutoff at
        # zero time, i.e., time of merger.
        ringing_crop = wave_p.start_time + wave_p.duration
        wave_p = wave_p.crop(0, ringing_crop)
        wave_c = wave_c.crop(0, ringing_crop)
        # Copy the unshifted waveform to save it as the template for matched-filtering
        self.mf_template = deepcopy(wave_p)
        # Shift merger time to 'time of explosion from detected sources'
        wave_p.start_time += gps_time_mod_year
        wave_c.start_time += gps_time_mod_year
        self.gw_signal_plus = wave_p
        self.gw_signal_cross = wave_c
