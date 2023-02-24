import unittest
from kilopop import mappings
from kilopop.kilonovae import bns_kilonova as kilonova


class test_mappings(unittest.TestCase):
    def test_equation_4(self):
        result = mappings.compute_equation_4(2.0, 1.5, 0.2, 0.15)
        self.assertAlmostEqual(result, 0.0057, None, 'Should be ', 0.0001)

    def test_equation_5(self):
        result = mappings.compute_equation_5(1.9, 1.8, 0.18, 0.15)
        self.assertAlmostEqual(result, 0.2310, None, 'Should be ', 0.0001)

    def test_equation_6(self):
        result = mappings.compute_equation_6(3.0, 3.5)
        self.assertAlmostEqual(result, 0.1397, None, 'Should be ', 0.0001)

    def test_equation_7(self):
        EOS_mass_to_rad = kilonova().__class__.EOS_mass_to_rad
        result = mappings.compute_equation_7(2.06, EOS_mass_to_rad)
        self.assertAlmostEqual(result, 3.601, None, 'Should be ', 0.001)

    def test_equation_8(self):
        result = mappings.compute_equation_8(0.2, 0.2)
        self.assertAlmostEqual(result, 0.0400, None, 'Should be ', 0.0001)

    def test_equation_9(self):
        result = mappings.compute_equation_9(0.002, 0.05)
        self.assertAlmostEqual(result, 0.0520, None, 'Should be ', 0.0001)

    def test_equation_10(self):
        result = mappings.compute_equation_10(0.5)
        self.assertAlmostEqual(result, 0.3363, None, 'Should be ', 0.0001)

    def test_gaussian_process_emulator(self):
        kn_inst = kilonova()
        gp = kn_inst.__class__.grey_opacity_emulator
        data = kn_inst.__class__.opacity_data
        grey_opacity = mappings.emulate_grey_opacity_from_kilonova_ejecta(0.05, 0.22, 0.35, gp, data)
        self.assertAlmostEqual(grey_opacity, 5.0, None, 'Should be ', 5.0)


if __name__ == '__main__':
    unittest.main()
