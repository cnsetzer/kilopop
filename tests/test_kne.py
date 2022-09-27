"""
 Module to test implementations of kne.py module
"""

import unittest
from bnspopkne.kne import Setzer2022_kilonova as kilonova
from bnspopkne.kne import Setzer2022_kilonova as population


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

    @unittest.skipIf
    def test_derived_properties_of_populations


if __name__ == '__main__':
    unittest.main()
