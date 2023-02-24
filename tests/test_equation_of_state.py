"""
 Module to test the functionality of handling the equation of state in equation_of_state.py
"""

import unittest
from pkg_resources import resource_filename
from kilopop import equation_of_state as eos
from kilopop.kilonovae import bns_kilonova as kilonova


class test_equation_of_state(unittest.TestCase):
    def setUp(self):
        self.EOS_table = eos.get_EOS_table(EOS_path=resource_filename('kilopop', "data/mr_sfho_full_right.csv"))
        self.radius_interpolator = eos.get_radius_interpolator_from_EOS(self.EOS_table)

    def test_radius_interpolator_from_EOS(self):
        test_radius = self.radius_interpolator(1.6)
        self.assertAlmostEqual(test_radius, 11.761, None, "Should be", 0.001)

    def test_max_EOS_mass(self):
        result = eos.get_max_EOS_mass(self.EOS_table)
        self.assertAlmostEqual(result, 2.056, None, "Should be", 0.001)

    def test_compactness_from_EOS(self):
        result = eos.compute_compactnesses_from_EOS(1.6, self.radius_interpolator)
        self.assertAlmostEqual(result, 0.201, None, "Placeholder", 0.001)

if __name__ == '__main__':
    unittest.main()
