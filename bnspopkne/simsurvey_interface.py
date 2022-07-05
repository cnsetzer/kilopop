"""Interface to use this kilonovae model with simsurvey for making observations."""
import numpy as np
import os
import sncosmo
import simsurvey
from astropy.cosmology import Planck18
from bnspopkne.kne import saee_bns_emgw_with_viewing_angle as saeev
from bnspopkne.population import Setzer2022_population as s2p


class saee_bns_emgw_with_viewing_angle_simsurvey(saeev):
    def __init__(self, **kwargs):
        pass

    def set(self, **kwargs):
        super().__init__(kwargs)

# Define the function that generates the lightcurve parameters
# (Note how the value for mags is not in the typical range for SNe Ia.
#  This will be fixed by passing new


def Setzer2022_population_simsurvey(
    redshifts,
    model,
    ra,
    dec,
    **kwargs,
):
    """
    """
    out =
    spop = s2p(num_samples=len(redshifts), **kwargs)
    for i in range(spop.num_params):
        out[getattr(spop, f"param{i}_name")] = getattr(spop, f"param{i}")

    out['EOS'] = np.array(list(repeat(kwargs['EOS'], len(redshifts))))
    out['EOS_path'] = np.array(list(repeat(kwargs['EOS_path'], len(redshifts))))
    out['kappa_grid_path'] = np.array(list(repeat(kwargs['kappa_grid_path'], len(redshifts))))
    out['gp_hyperparameter_file'] = np.array(list(repeat(kwargs['gp_hyperparameter_file'], len(redshifts))))
    out['mapping_type'] = np.array(list(repeat(kwargs['mapping_type'], len(redshifts))))
    out['cosmo'] = np.array(list(repeat(Planck18, len(redshifts))))
    out['z'] = redshifts
    out['ra'] = ra
    out['dec'] = dec
    return out


def simulate_simsurvey_population(zmin, zmax, rate, plan):
    transientprop = {
        "lcmodel": saee_bns_emgw_with_viewing_angle_simsurvey,
        "lcsimul_func": Setzer2022_population_simsurvey,
        "lcsimul_prop": {'EOS': 'sfho',
                         'EOS_path': None,
                         'gp_hyperparameter_file': None,
                         'kappa_grid_path': None,
                         'mapping_type': "coughlin",
                         }}

    tr = simsurvey.get_transient_generator(
        (zmin, zmax),
        ratefunc=lambda z: rate,
        ra_range=(0.0, 360.0),
        dec_range=(-90.0, 15.0),
        mjd_range=(58178, 58543),
        transientprop=transientprop,
    )
    survey = simsurvey.SimulSurvey(generator=tr, plan=plan)
    lcs = survey.get_lightcurves(progress_bar=True)
    # The lightcurves can further be saved as a pickle file
    lcs.save("lcs.pkl")
