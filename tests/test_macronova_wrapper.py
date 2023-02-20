import unittest
import pickle
import numpy as np
from pkg_resources import resource_filename
from bnspopkne import macronovae_wrapper as mw


class test_macronova_engine(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(resource_filename('bnspopkne', "data/KNE_parameters.pkl"), "rb") as input:
            cls._kne_parameters = pickle.load(input)
        with open(resource_filename('bnspopkne', "data/phase.pkl"), "rb") as input:
            cls._known_phase = pickle.load(input)
        with open(resource_filename('bnspopkne', "data/wave.pkl"), "rb") as input:
            cls._known_wave = pickle.load(input)
        with open(resource_filename('bnspopkne', "data/flux.pkl"), "rb") as input:
            cls._known_flux = pickle.load(input)

    def test_planck_function(self):
        # checking peak value is reproduced
        peak_10k = 4095673082270000.0  # erg/s /sr /cm^2 /Ang
        lambda_cm = np.array([2.898 * 1e-5])  # wavelength of 3000 Ang in cm
        temperature = np.array([10000])  # K
        result = mw.compute_planck_function(lambda_cm, temperature)
        # checking if within .01% of the known value
        self.assertAlmostEqual(result, peak_10k, None, "", 0.0001 * peak_10k)

    def test_create_seds(self):
        """Test the SED creation functionality.
        Inadvertenty this also tests create_sed_timeseries.
        """
        # load known solution
        phase, wave, flux = mw.create_saee_seds(self._kne_parameters)
        np.testing.assert_allclose(phase, self._known_phase, 0.01, err_msg='Should be equal:', verbose=True)
        np.testing.assert_allclose(wave, self._known_wave, 0.01, err_msg='Should be equal:', verbose=True)
        np.testing.assert_allclose(flux, self._known_flux, 0.01, err_msg='Should be equal:', verbose=True)


if __name__ == '__main__':
    unittest.main()
