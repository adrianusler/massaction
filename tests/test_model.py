from massaction.model import ChemModel
from massaction.species import Species

from unittest import TestCase


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
