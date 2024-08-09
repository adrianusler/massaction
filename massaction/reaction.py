"""Module for the definition of chemical reactions."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from massaction.model import SpeciesLike

import numpy as np


class Reaction:
    """Class for modeling chemical reactions.

    :param reactants: LinCombSpecies object representing the reactants.
    :param products: LinCombSpecies object representing the products.
    """

    def __init__(self, reactants: SpeciesLike, products: SpeciesLike) -> None:
        """Instantiate Reaction object."""
        from massaction.species import ensure_lincomb

        self.reactants = ensure_lincomb(reactants)
        self.products = ensure_lincomb(products)
        self.lincomb = self.products - self.reactants

    def eval(self, ln_concentrations: np.ndarray, ln_eqconst: float) -> float:
        """Evaluate the (logarithmic) mass-action law of the reaction.

        :param ln_concentrations: Array of logarithmic concentrations.
        :param ln_eqconst: logarithmic equilibrium constant.
        """
        factors = np.array([fs.factor for fs in self.lincomb.factor_species_list])
        species_ids = [fs.species_id for fs in self.lincomb.factor_species_list]
        return np.sum(factors * ln_concentrations[species_ids]) - ln_eqconst

    def print(self) -> None:
        """Print the specifics of the Reaction object."""
        reactants_str = self.reactants.print(return_str=True)
        products_str = self.products.print(return_str=True)
        print(f"{reactants_str} >> {products_str}")  # noqa: T201
