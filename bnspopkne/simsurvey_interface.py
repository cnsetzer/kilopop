"""Interface to use this kilonovae model with simsurvey for making observations."""
import numpy as np
import os
import sncosmo
import simsurvey
from astropy.cosmology import Planck15


# Define the function that generates the lightcurve parameters
# (Note how the value for mags is not in the typical range for SNe Ia.
#  This will be fixed by passing new value
def random_parameters(
    redshifts,
    model,
    mag=(-18.0, 1.0),
    r_v=2.0,
    ebv_rate=0.11,
    alpha=1.3,
    cosmo=Planck15,
    **kwargs
):
    """
    """
    out = {}

    amp = []
    for z in redshifts:
        mabs = np.random.normal(mag[0], mag[1])
        model.set(z=z)
        model.set_source_peakabsmag(mabs, "bessellb", "vega", cosmo=cosmo)
        amp.append(model.get("amplitude"))

    out["amplitude"] = np.array(amp)
    out["hostr_v"] = r_v * np.ones(len(redshifts))
    out["hostebv"] = np.random.exponential(ebv_rate, len(redshifts))

    out["s"] = np.random.normal(1.0, 0.1, len(redshifts))
    out["amplitude"] *= 10 ** (0.4 * alpha * (out["s"] - 1))

    return out


transientprop = {
    "lcmodel": model,
    "lcsimul_func": random_parameters,
    "lcsimul_prop": {"mag": (-19.3, 0.1)},
}

tr = simsurvey.get_transient_generator(
    (0.0, 0.05),
    ratefunc=lambda z: 3e-5,
    ra_range=(0, 360),
    dec_range=(-30, 90),
    mjd_range=(58178, 58543),
    transientprop=transientprop,
)

plan = None  # Place holder to get simsurvey-like plan from OpSim


survey = simsurvey.SimulSurvey(generator=tr, plan=plan)

lcs = survey.get_lightcurves(progress_bar=True)

# The lightcurves can further be saved as a pickle file
lcs.save("lcs.pkl")
