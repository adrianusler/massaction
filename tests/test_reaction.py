from massaction.model import ChemModel
from massaction.species import LinCombSpecies
from massaction.reaction import Reaction

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

    def test_reaction(self):
        h2o, h2, o2 = self.h2o, self.h2, self.o2
        reaction = self.reaction

        assert isinstance(reaction, Reaction)
        assert reaction.reactants == h2 + o2
        assert reaction.products == LinCombSpecies([2 * h2o])

        ln_conc = np.array([10.6, 11.2, 9.3])
        ln_eqconst = -7.0
        ln_massaction = 2 * ln_conc[0] - ln_conc[1] - ln_conc[2] - ln_eqconst
        assert np.isclose(reaction.eval(ln_conc, ln_eqconst), ln_massaction)

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def test_print_reaction(self):
        reaction = self.reaction
        reaction.print()
        captured_text = self.capsys.readouterr().out
        assert captured_text == "+1.0*species_1 +1.0*species_2  >> +2.0*species_0 \n"
