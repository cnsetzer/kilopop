import unittest
from bnspopkne.kne import Setzer2022_kilonova as kilonova
from bnspopkne.kne import Setzer2022_kilonova as population


class test_kne(unittest.TestCase):
    def test_kilonova(self):
        kilo_instance = kilonova()

    def test_kilonova_population(self):
        pop_instance = population()


if __name__ == '__main__':
    unittest.main()
