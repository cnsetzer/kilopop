"""
 Module to test implementations of kne.py module
"""

import unittest
import numpy as np
import sncosmo
from kilopop.kilonovae import bns_kilonova as kilonova
from kilopop.kilonovae import bns_kilonovae_population_distribution as population


class test_kne(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._paper_ranges = {'param1':(1.0, 2.06),'param2':(1.0,2.06),'param3':(0.1,0.3),'param4':(0.1,0.3),'param5':(0.0,2*np.pi),'param6':(0.15,0.4),'param7':(0.001,0.08),'param8':(0.17,0.38),'param9':(0.01,500.0),'param10':(0.0001,0.08),'param11':(0.001,0.08),'param12':(0.1,0.4)}
        cls._paper_derived = {'peak_time':(0.1, 13.0), 'peak_absmag_lssti':(-20.3,-10.5),'one_mag_peak_time_lssti':(0.0,15.0)}
        _ = sncosmo.get_bandpass('lssti')

    def test_kilonova(self):
        kilo_instance = kilonova()
        self.assertIsInstance(kilo_instance, kilonova)
        for key, value in self._paper_ranges.items():
            inst_val = getattr(kilo_instance, key)
            self.assertTrue((value[0] <= inst_val) & (inst_val <= value[1]), msg=f"The value for {key} is {inst_val}, which is not in the range of {value[0]} to {value[1]}")


    def test_kilonova_population(self):
        pop_instance = population()
        for key, value in self._paper_ranges.items():
            inst_val = getattr(pop_instance, key)
            self.assertTrue(all((value[0] <= inst_val) & (inst_val <= value[1])), msg=f"The value for {key} is {inst_val}, which is not in the range of {value[0]} to {value[1]}")

    def test_derived_properties_of_populations(self):
        pop_instance = population(population_size=1000, only_draw_parameters=False)
        for key, value in self._paper_ranges.items():
            inst_val = getattr(pop_instance, key)
            self.assertTrue(all((value[0] <= inst_val) & (inst_val <= value[1])), msg=f"The value for {key} is {inst_val}, which is not in the range of {value[0]} to {value[1]}")
        for key, value in self._paper_derived.items():
            inst_val = getattr(pop_instance, key)
            self.assertTrue(all((value[0] <= inst_val) & (inst_val <= value[1])), msg=f"The value for {key} is {inst_val}, which is not in the range of {value[0]} to {value[1]}")


if __name__ == '__main__':
    unittest.main()
