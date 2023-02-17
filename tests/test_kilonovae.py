"""
 Module to test implementations of kne.py module
"""

import unittest
from bnspopkne.kilonovae import Setzer2022_kilonova as kilonova
from bnspopkne.kilonovae import Setzer2022_population_parameter_distribution as population


class test_kne(unittest.TestCase):
    def test_kilonova(self):
        kilo_instance = kilonova()
        self.assertIsInstance(kilo_instance, kilonova)

    def test_kilonova_population(self):
        paper_ranges = {}
        paper_properties = {}
        pop_instance = population()
        for key, value in paper_ranges.items():
            getattr(pop_instance, key)

    def test_derived_properties_of_populations(self):
        pop_instance = population(population_size=1000, only_draw_parameters=False)
        paper_ranges = {}


if __name__ == '__main__':
    unittest.main()
