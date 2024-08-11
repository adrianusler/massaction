from massaction.model import ChemModel

from unittest import TestCase
import numpy as np

import pytest


class TestReaction(TestCase):
    """Test the Reaction class."""

    def setUp(self):
        self.model3 = ChemModel(3)
        h2o, h2, o2 = self.model3.get_all_species()
        self.h2o, self.h2, self.o2 = h2o, h2, o2
        self.reaction = h2 + o2 >> 2 * h2o

    def test_reaction_setups(self):
        def reaction_assertions(
            reaction, reactants, reactants_factors, products, products_factors
        ):
            lnc = np.array([-0.3, 1.5, 2.7])
            lnK = 4.2
            mal_eval = -lnK
            for i, reactant in enumerate(reactants):
                fs = reaction.reactants.factor_species_list[i]
                assert fs.species is reactant
                assert fs.factor == reactants_factors[i]
                mal_eval -= reactants_factors[i] * lnc[reactant.species_id]
            for i, product in enumerate(products):
                fs = reaction.products.factor_species_list[i]
                assert fs.species is product
                assert fs.factor == products_factors[i]
                mal_eval += products_factors[i] * lnc[product.species_id]
            assert np.isclose(reaction.eval(lnc, lnK), mal_eval)

        h2o, h2, o2 = self.h2o, self.h2, self.o2
        # test Species.__rshift__
        reaction = h2 >> 1.9 * o2
        reaction_assertions(reaction, [h2], [1.0], [o2], [1.9])
        # test FactorSpecies.__rshift__
        reaction = 0.4 * h2 >> o2
        reaction_assertions(reaction, [h2], [0.4], [o2], [1.0])
        # test LinCombSpecies.__rshift__
        reaction = (h2 + o2) >> 2 * h2o
        reaction_assertions(reaction, [h2, o2], [1.0, 1.0], [h2o], [2.0])

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def test_print_reaction(self):
        reaction = self.reaction
        reaction.print()
        captured_text = self.capsys.readouterr().out
        assert captured_text == "+1.0*species_1 +1.0*species_2  >> +2.0*species_0 \n"
