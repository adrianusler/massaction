from massaction.model import ChemModel
from massaction.species import LinCombSpecies
from massaction.reaction import Reaction

from unittest import TestCase
import numpy as np


class TestReaction(TestCase):
    """Test the Reaction class."""

    def setUp(self):
        self.model3 = ChemModel(3)

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
