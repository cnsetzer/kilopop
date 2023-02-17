import unittest
import numpy as np
from bnspopkne import population_priors


class test_population(unittest.TestCase):
    def test_disk_unbinding_efficiency(self):
        result = population_priors.draw_disk_unbinding_efficiency(output_shape=10)
        for res in result:
            self.assertTrue(0.1 <= res <= 0.4)

    def test_viewing_angle(self):
        result = population_priors.draw_viewing_angle(output_shape=10)
        for res in result:
            self.assertTrue(0.0 <= res <= np.pi/2.0)

    def test_mass_from_EOS_bounds(self):
        result = population_priors.draw_mass_from_EOS_bounds(2.05, output_shape=10)
        for res in result:
            self.assertTrue(1.0 <= res <= 2.05)

    def test_mass_from_EOS_bounds_with_mass_ratio_cut(self):
        result1, result2 = population_priors.draw_masses_from_EOS_bounds_with_mass_ratio_cut(2.05, mass_ratio_cut=2.0/3.0, output_shape=10)
        for res in result1:
            self.assertTrue(1.0 <= res <= 2.05)
        for i, res2 in enumerate(result2):
            print(max(result1[i]*2.0/3.0, 1.0), res2, result1[i])
            self.assertTrue(max(result1[i]*2.0/3.0, 1.0) <= res2 <= result1[i])


if __name__ == '__main__':
    unittest.main()
