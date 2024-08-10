from massaction.model import ChemModel
from massaction.species import FactorSpecies, LinCombSpecies

from unittest import TestCase


class TestSpecies(TestCase):
    """Test the SpeciesLike classes Species, FactorSpecies, and LinCombSpecies."""

    def setUp(self):
        self.model3 = ChemModel(3)

    def test_lincomb(self):
        h2o, h2, o2 = self.model3.get_all_species()
        h2o_2 = 2 * h2o
        assert isinstance(h2o_2, FactorSpecies)
        assert h2o_2.factor == 2.0

        lincomb1 = 3.4 * h2o - 0.2 * h2 - o2
        lincomb2 = o2 - 2.0 * h2o + h2
        assert isinstance(lincomb1, LinCombSpecies)
        assert isinstance(lincomb2, LinCombSpecies)
        assert len(lincomb1.factor_species_list) == 3
        assert len(lincomb2.factor_species_list) == 3

        factors1 = [3.4, -0.2, -1.0]
        factors2 = [-2.0, 1.0, 1.0]
        for fs in lincomb1.factor_species_list:
            assert isinstance(fs, FactorSpecies)
            assert fs.factor == factors1[fs.species_id]
        for fs in lincomb2.factor_species_list:
            print(fs.species_id, fs.factor)
            assert isinstance(fs, FactorSpecies)
            assert fs.factor == factors2[fs.species_id]
