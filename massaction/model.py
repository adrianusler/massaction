"""Module for defining mass-action laws with constraints."""

from __future__ import annotations

import numpy as np
from scipy.optimize import root

from massaction.constraint import Constraint
from massaction.reaction import Reaction
from massaction.species import Species


class ChemModel:
    """Class for modeling chemical species. Manages species and reactions."""

    def __init__(self, num_species: int) -> None:
        """Instantiate ChemModel object."""
        self.num_species = num_species
        self.list_of_species = [Species(i, self) for i in range(num_species)]
        self.list_of_reactions = []
        self.list_of_constraints = []

    def get_all_species(self) -> list:
        """Return list of all species."""
        return self.list_of_species

    def solve(
        self,
        reactions: list[Reaction],
        ln_eqconst: list[float],
        constraints: list[Constraint],
    ) -> np.ndarray:
        """Solve the system of reactions and constraints.

        :param reactions: List of Reaction objects.
        :param ln_eqconst: Array of logarithmic equilibrium constants.
        :param constraints: List of Constraint objects.
        """
        num_reactions = len(reactions)
        num_constraints = len(constraints)
        num_species = self.num_species
        num_eq = num_reactions + num_constraints
        if num_eq != num_species:
            msg = "Number of equations must be equal to the number of species."
            raise ValueError(msg)
        if len(ln_eqconst) != num_reactions:
            msg = "Number of equilibrium constants must be equal to the number of reactions."
            raise ValueError(msg)

        def eval_system_of_equations(ln_concentrations: np.ndarray) -> np.ndarray:
            ln_equation = np.zeros(num_eq)
            for i in range(num_reactions):
                ln_equation[i] = reactions[i].eval(ln_concentrations, ln_eqconst[i])
            for i in range(num_constraints):
                ln_equation[num_reactions + i] = constraints[i].eval(ln_concentrations)
            return ln_equation

        first_guess = np.zeros(num_species)
        result = root(eval_system_of_equations, x0=first_guess)

        return result.x
