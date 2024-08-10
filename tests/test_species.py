from massaction.model import ChemModel
from massaction.species import FactorSpecies, LinCombSpecies

from unittest import TestCase


class TestSpecies(TestCase):
    """Test the SpeciesLike classes Species, FactorSpecies, and LinCombSpecies."""

    def setUp(self):
        self.model3 = ChemModel(3)
        self.h2o, self.h2, self.o2 = self.model3.get_all_species()

    def test_factor_species(self):
        def fs_assertions(fs, factor, species_obj):
            assert isinstance(fs, FactorSpecies)
            assert fs.factor == factor
            assert fs.species == species_obj

        h2o, h2, o2 = self.h2o, self.h2, self.o2
        # test Species.__mul__
        fs_assertions(2.0 * h2o, 2.0, h2o)
        # test Species.__neg__
        fs_assertions(-o2, -1.0, o2)
        # test FactorSpecies.__neg__
        fs_assertions(-(6.0 * h2), -6.0, h2)
        # test FactorSpecies.__mul__
        fs_assertions(4.0 * (-3 * h2o), -12.0, h2o)

    def test_lincomb(self):
        def lc_assertions(lc, factors_list, species_list):
            assert isinstance(lc, LinCombSpecies)
            assert len(factors_list) == len(
                species_list
            )  # just internal consistency check
            assert len(lc.factor_species_list) == len(factors_list)
            lc_factors = [fs.factor for fs in lc.factor_species_list]
            lc_species = [fs.species for fs in lc.factor_species_list]
            assert lc_factors == factors_list
            for species1, species2 in zip(lc_species, species_list):
                assert species1.species_id == species2.species_id

        h2o, h2, o2 = self.h2o, self.h2, self.o2

        # test Species.__add__
        lc_assertions(h2o + o2, [1.0, 1.0], [h2o, o2])
        # test Species.__sub__
        lc_assertions(h2 - h2o, [1.0, -1.0], [h2, h2o])
        # test FactorSpecies.__add__
        lc_assertions(1.3 * h2o + o2, [1.3, 1.0], [h2o, o2])
        # test FactorSpecies.__sub__
        lc_assertions(1.3 * h2o - 7.0 * o2, [1.3, -7.0], [h2o, o2])
        lc_assertions(1.3 * h2o - (-3.1) * o2, [1.3, 3.1], [h2o, o2])
        # test LinComb.__add__
        lc_assertions((-3.0 * h2 + o2) + 0.6 * h2o, [-3.0, 1.0, 0.6], [h2, o2, h2o])
        # test LinComb.__sub__
        lc_assertions(
            (4.0 * h2o - 11.0 * o2) - 21 * h2, [4.0, -11.0, -21.0], [h2o, o2, h2]
        )
