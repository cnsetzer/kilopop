import numpy as np
import astropy.units as uu
import sncosmo
import redback
from bnspopkne.kilonovae import Setzer2022_kilonova


class redback_wrapper(Object):

    def __init__(time, redshift, band, **kwargs):
        kilonova = Setzer2022_kilonova(redshift=redshift, kwargs)
        ff = kilonova.model.bandflux(time, band)
        ff = ff * uu.erg/uu.s/(uu.cm*uu.cm)

        if kwargs['output_format'] == 'flux':
            return ff.to(uu.mJy).value
        elif kwargs['output_format'] == 'magnitude':
            return ff.to(uu.ABmag).value

        # convert to source frame time and frequency
        time = time * day_to_s
        frequency, time = calc_kcorrected_properties(frequency=frequency, redshift=redshift, time=time)
        dl = cosmo.luminosity_distance(redshift).cgs.value
