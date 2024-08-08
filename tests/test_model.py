from lib.defchem import (
    ChemModel,
    Species,
    FactorSpecies,
    LinCombSpecies,
    Constraint,
    Reaction,
)
from unittest import TestCase
import numpy as np


class TestChemModel(TestCase):
    """Test the ChemModel class."""

    def setUp(self):
        self.model1 = ChemModel(1)
        self.model2 = ChemModel(2)
        self.model3 = ChemModel(3)

    def test_get_all_species(self):
        """Test the get_all_species method."""
        list_species = self.model1.get_all_species()
        assert isinstance(list_species, list)
        assert len(list_species) == 1
        assert all(isinstance(s, Species) for s in list_species)
        list_species = self.model2.get_all_species()
        assert isinstance(list_species, list)
        assert len(list_species) == 2
        assert all(isinstance(s, Species) for s in list_species)
        list_species = self.model3.get_all_species()
        assert isinstance(list_species, list)
        assert len(list_species) == 3
        assert all(isinstance(s, Species) for s in list_species)

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

    def test_constraint(self):
        h2o, h2, o2 = self.model3.get_all_species()

        lincomb = 3.4 * h2o - 0.2 * h2 - o2
        constraint = lincomb == -1.3
        assert isinstance(constraint, Constraint)
        assert constraint.lincomb == lincomb
        assert constraint.value == -1.3
        lnconc = np.array([23.0, -1.0, 4.3])
        conc = np.exp(lnconc)
        constraint_eval_pos = 3.4 * conc[0] + 1.3
        constraint_eval_neg = 0.2 * conc[1] + conc[2]
        assert np.isclose(
            constraint.eval(lnconc),
            np.log(constraint_eval_pos) - np.log(constraint_eval_neg),
        )

        constraint = lincomb == 10.7
        assert isinstance(constraint, Constraint)
        assert constraint.lincomb == lincomb
        assert constraint.value == 10.7
        lnconc = np.array([3.14, 3.0, 2.9])
        conc = np.exp(lnconc)
        constraint_eval_pos = 3.4 * conc[0]
        constraint_eval_neg = 0.2 * conc[1] + conc[2] + 10.7
        assert np.isclose(
            constraint.eval(lnconc),
            np.log(constraint_eval_pos) - np.log(constraint_eval_neg),
        )

    def test_reaction(self):
        h2o, h2, o2 = self.model3.get_all_species()
        reaction = h2 + o2 >> 2 * h2o
        assert isinstance(reaction, Reaction)
        assert reaction.reactants == h2 + o2
        assert reaction.products == LinCombSpecies([2 * h2o])

        ln_conc = np.array([10.6, 11.2, 9.3])
        ln_eqconst = -7.0
        ln_massaction = 2 * ln_conc[0] - ln_conc[1] - ln_conc[2] - ln_eqconst
        assert np.isclose(reaction.eval(ln_conc, ln_eqconst), ln_massaction)
