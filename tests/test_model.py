from massaction.model import ChemModel, get_num_sweep
from massaction.species import Species

from unittest import TestCase

import pytest
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

    def test_get_num_sweep(self):
        h2o, h2, o2 = self.model3.get_all_species()
        c1 = h2o == 1.0
        c2 = h2 == 4.0
        c3 = o2 == [1.0, 2.0, 3.0, 4.0, 5.0]
        assert get_num_sweep([c1, c2, c3]) == 5
        c2 = h2 == [10.0, 11.0, 12.0]
        with pytest.raises(ValueError):
            get_num_sweep([c1, c2, c3])

    def test_solve(self):
        h2o, h2, o2 = self.model3.get_all_species()
        reaction = h2 + o2 >> 2 * h2o
        ln_eqconst = 3.0
        constraint1 = h2o + h2 == 1.0
        constraint2 = h2o + 2 * o2 == 7.0
        with pytest.raises(ValueError):
            self.model3.solve([reaction], [ln_eqconst], [constraint1])
        with pytest.raises(ValueError):
            self.model3.solve(
                [reaction], [ln_eqconst, ln_eqconst], [constraint1, constraint2]
            )
        lnc = self.model3.solve([reaction], [ln_eqconst], [constraint1, constraint2])
        assert np.isclose(2 * lnc[0] - lnc[1] - lnc[2], ln_eqconst)
        c = np.exp(lnc)
        assert np.isclose(c[0] + c[1], 1.0)
        assert np.isclose(c[0] + 2 * c[2], 7.0)
