"""Module for defining mass-action laws with constraints."""

from __future__ import annotations


import numpy as np
from scipy.optimize import root

from massaction.constraint import Constraint, ConstraintSweep, ConstraintLike
from massaction.reaction import Reaction
from massaction.species import Species, LinCombSpecies


class ChemModel:
    """Class for modeling chemical species. Manages species and reactions."""

    def __init__(self, num_species: int) -> None:
        """Instantiate ChemModel object."""
        self.num_species = num_species
        self.list_of_species = [Species(i, self) for i in range(num_species)]

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
        # no. of equations (including reservoir constraints)
        num_eq_total = num_reactions + num_constraints
        if num_eq_total != num_species:
            msg = "Number of equations must be equal to the number of species."
            raise ValueError(msg)
        if len(ln_eqconst) != num_reactions:
            msg = "Number of equilibrium constants must be equal to the number of reactions."
            raise ValueError(msg)

        # Check for reservoir constraints (constraints with only one species)
        check_reservoirs = [cstr.check_reservoir() for cstr in constraints]
        reservoir_species_ids = [
            check_res[1] for check_res in check_reservoirs if check_res[0]
        ]
        non_reservoir_species_ids = [
            i for i in range(num_species) if i not in reservoir_species_ids
        ]
        reservoir_constraint_ids = [
            i for i, check_res in enumerate(check_reservoirs) if check_res[0]
        ]
        non_reservoir_constraint_ids = [
            i for i, check_res in enumerate(check_reservoirs) if not check_res[0]
        ]

        num_reservoirs = len(reservoir_species_ids)
        num_eq = num_eq_total - num_reservoirs

        def update_reservoir_values():
            """Update the reservoir values."""
            reservoir_values[:] = [
                constraints[cstr_id].reservoir_value
                for cstr_id in reservoir_constraint_ids
            ]

        reservoir_values = [0.0] * num_reservoirs
        update_reservoir_values()

        def ln_concentrations_with_reservoirs(
            ln_concentrations: np.ndarray,
        ) -> np.ndarray:
            """Add reservoir concentrations to the logarithmic concentrations."""
            ln_all_concentrations = np.zeros(num_species)
            ln_all_concentrations[reservoir_species_ids] = np.log(reservoir_values)
            ln_all_concentrations[non_reservoir_species_ids] = ln_concentrations
            return ln_all_concentrations

        def eval_system_of_equations(ln_concentrations: np.ndarray) -> np.ndarray:
            """Evaluate the system of equations."""
            ln_all_concentrations = ln_concentrations_with_reservoirs(ln_concentrations)
            ln_equation = np.zeros(num_eq)

            for i in range(num_reactions):
                ln_equation[i] = reactions[i].eval(ln_all_concentrations, ln_eqconst[i])
            for i in range(num_constraints - num_reservoirs):
                constraint_id = non_reservoir_constraint_ids[i]
                ln_equation[num_reactions + i] = constraints[constraint_id].eval(
                    ln_all_concentrations
                )
            return ln_equation

        num_sweep = get_num_sweep(constraints)
        _num_sweep = 1 if num_sweep == -1 else num_sweep
        results_list = []
        for i in range(_num_sweep):
            [cstr.set_current_id(i) for cstr in constraints]
            update_reservoir_values()
            first_guess = np.zeros(num_species - num_reservoirs)
            result = root(eval_system_of_equations, x0=first_guess)
            results_list += [ln_concentrations_with_reservoirs(result.x)]
        if num_sweep == -1:
            return results_list[0]
        return np.array(results_list)


def get_num_sweep(constraints: list[ConstraintLike]) -> int:
    """Return the number of values in the parameter sweep. If inconsistent, raise error."""
    num_sweep = -1
    for cstr in constraints:
        if not isinstance(cstr, ConstraintSweep):
            continue
        if num_sweep == cstr.num_values:
            continue
        if num_sweep == -1:
            num_sweep = cstr.num_values
            continue
        msg = "All constraint sweeps must have equal number of values"
        raise ValueError(msg)
    return num_sweep


nil = LinCombSpecies([])
