"""Module for the definition of chemical constraints."""

from __future__ import annotations

from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from massaction.model import SpeciesLike

import numpy as np


class Constraint:
    """Class for modeling constraints on chemical species."""

    def __init__(self, lincomb: SpeciesLike, value: float) -> None:
        """Instantiate Constraint object.

        :param lincomb: LinCombSpecies object.
        :param value: Value of the constraint.
        """
        from massaction.species import ensure_lincomb

        self.lincomb = ensure_lincomb(lincomb)
        self.value = value

    def eval(self, ln_concentrations: np.ndarray) -> float:
        """Evaluate the constraint.

        :param ln_concentrations: Array of logarithmic concentrations.
        """
        concentrations = np.exp(ln_concentrations)
        factors = np.array([fs.factor for fs in self.lincomb.factor_species_list])
        species_ids = [fs.species_id for fs in self.lincomb.factor_species_list]
        mask_pos = factors > 0.0
        mask_neg = factors < 0.0
        val_pos = (
            np.sum(factors[mask_pos] * concentrations[species_ids][mask_pos])
            if mask_pos.any()
            else 0.0
        )
        val_neg = (
            -np.sum(factors[mask_neg] * concentrations[species_ids][mask_neg])
            if mask_neg.any()
            else 0.0
        )
        if self.value > 0:
            return np.log(val_pos) - np.log(val_neg + self.value)
        return np.log(val_pos + (-self.value)) - np.log(val_neg)

    def print(self) -> None:
        """Print the specifics of the Constraint object."""
        mystr = self.lincomb.print(return_str=True)
        print(f"{mystr} == {self.value}")  # noqa: T201

    @property
    def num_values(self) -> int:
        return 1

    def set_current_id(self, current_id: int) -> None:
        pass


class ConstraintSweep(Constraint):
    def __init__(self, lincomb: SpeciesLike, values: list | np.ndarray):
        self.values = values
        self.current_id = 0
        super().__init__(lincomb, self.values[self.current_id])

    @property
    def num_values(self) -> int:
        return len(self.values)

    def set_current_id(self, current_id: int) -> None:
        self.current_id = current_id
        self.value = self.values[self.current_id]


ConstraintLike = Union[Constraint, ConstraintSweep]
