import unittest
from bnspopkne import macronovae_wrapper as mw


class test_macronova_engine(unittest.TestCase):
    def setUp(self):

    def test_planck_function():
        lamdba_cm = np.array(3.0*1e-5)  # wavelength of 3000 Ang in cm
        temperature = np.array(10000)  # K
        result = mw.compute_planck_function(lambda_cm, temperature)
        self.assertAlmostEqual

    def test_sed_timeseries():


    def test_create_seds():



if __name__ == '__main__':
    unittest.main()
