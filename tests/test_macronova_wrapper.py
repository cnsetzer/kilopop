import unittest
import pickle
from pkg_resources import resource_filename
from bnspopkne import macronovae_wrapper as mw


class test_macronova_engine(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(resource_filename('bnspopkne', "data/KNE_parameters.pkl"), "rb") as input:
            cls._kne_parameters = pickle.load(input)

    def test_planck_function(self):
        # checking peak value is reproduced
        peak_10k = 4095673082270000.0  # erg/s /sr /cm^2 /Ang
        lamdba_cm = np.array([2.898*1e-5])  # wavelength of 3000 Ang in cm
        temperature = np.array([10000])  # K
        result = mw.compute_planck_function(lambda_cm, temperature)
        # checking if within .01% of the known value
        self.assertAlmostEqual(result, peak_10k, None, "", 0.0001*peak_10k)

    def test_create_seds(self, cls):
        """Test the SED creation functionality.
        Inadvertenty this also tests create_sed_timeseries.
        """
        # load known solution
        phase, wave, flux = mw.create_saee_seds(cls._kne_parameters)
        self.assertAlmostEqual(phase, known_phase, None, "", )
        self.assertAlmostEqual(wave, known_wave, None, "", )
        self.assertAlmostEqual(flux, known_flux, None, "", )


if __name__ == '__main__':
    unittest.main()
