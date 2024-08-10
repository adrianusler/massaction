"""Module for modelling chemical species."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from massaction.model import ChemModel
from numpy import ndarray

from massaction.constraint import Constraint, ConstraintSweep, ConstraintLike
from massaction.reaction import Reaction

from typing import Union


class Species:
    """Class for modeling chemical species."""

    def __init__(self, species_id: int, parent_model: ChemModel) -> None:
        """Instantiate Species object.

        :param species_id: Species ID number.
        :param parent_model: Parent ChemModel object.
        """
        self.species_id = species_id
        self.parent_model = parent_model
        self.ln_concentration = 0.0

    def __neg__(self) -> FactorSpecies:
        """Return FactorSpecies object representing the negation of a Species object."""
        return FactorSpecies(self, -1.0)

    def __mul__(self, factor: float) -> FactorSpecies:
        """Return FactorSpecies object representing the product of a species and a numerical factor."""
        return FactorSpecies(self, float(factor))

    def __rmul__(self, factor: float) -> FactorSpecies:
        """Return FactorSpecies object representing the product of a species and a numerical factor."""
        return FactorSpecies(self, float(factor))

    def __add__(self, other: SpeciesLike) -> LinCombSpecies:
        """Return LinCombSpecies object representing the sum of two Species objects."""
        return ensure_lincomb(self) + ensure_lincomb(other)

    def __sub__(self, other: SpeciesLike) -> LinCombSpecies:
        """Return LinCombSpecies object representing the difference of two Species objects."""
        return ensure_lincomb(self) - ensure_lincomb(other)

    def __eq__(self, value: Union[float, list[float], ndarray]) -> ConstraintLike:
        """Return Constraint object representing the equality of a Species object and a LinCombSpecies object."""
        if isinstance(value, Union[list, ndarray]):
            return ConstraintSweep(self, value)
        return Constraint(self, value)

    def __rshift__(self, other: SpeciesLike) -> Reaction:
        """Return Reaction object representing a chemical reaction."""
        return Reaction(self, other)


class FactorSpecies:
    """Class for modeling the union of a chemical species and a numerical pre-factor."""

    def __init__(self, species: Species, factor: float) -> None:
        """Instantiate FactorSpecies object.

        :param species: Species object.
        :param pre_factor: Numerical pre-factor.
        """
        self.species = species
        self.factor = factor

    @property
    def species_id(self) -> int:
        """Return species ID number."""
        return self.species.species_id

    def __neg__(self) -> FactorSpecies:
        """Return FactorSpecies object representing the negation of a FactorSpecies object."""
        return FactorSpecies(self.species, -1.0 * self.factor)

    def __mul__(self, factor: float) -> FactorSpecies:
        """Return FactorSpecies object representing the product of a FactorSpecies and a numerical factor."""
        return FactorSpecies(self.species, self.factor * factor)

    def __rmul__(self, factor: float) -> FactorSpecies:
        """Return FactorSpecies object representing the product of a FactorSpecies and a numerical factor."""
        return FactorSpecies(self.species, self.factor * factor)

    def __add__(self, other: SpeciesLike) -> LinCombSpecies:
        """Return LinCombSpecies object representing the sum of two FactorSpecies objects."""
        return ensure_lincomb(self) + ensure_lincomb(other)

    def __sub__(self, other: SpeciesLike) -> LinCombSpecies:
        """Return LinCombSpecies object representing the difference of two FactorSpecies objects."""
        return ensure_lincomb(self) - ensure_lincomb(other)

    def __eq__(self, value: Union[float, list[float], ndarray]) -> ConstraintLike:
        """Return Constraint object representing the equality of a FactorSpecies object and a numerical value."""
        if isinstance(value, Union[list, ndarray]):
            return ConstraintSweep(self, value)
        return Constraint(self, value)

    def __rshift__(self, other: SpeciesLike) -> Reaction:
        """Return Reaction object representing a chemical reaction."""
        return Reaction(self, other)


class LinCombSpecies:
    """Class for modeling the linear combination of chemical species."""

    def __init__(self, factor_species_list: list) -> None:
        """Instantiate LinCombSpecies object.

        :param factor_species_list: List of FactorSpecies objects.
        """
        self.factor_species_list = factor_species_list

    def __neg__(self) -> LinCombSpecies:
        """Return LinCombSpecies object representing the negation of a FactorSpecies object."""
        _factor_species_list = [-fs for fs in self.factor_species_list]
        return LinCombSpecies(_factor_species_list)

    def __add__(self, other: SpeciesLike) -> LinCombSpecies:
        """Return LinCombSpecies object representing the sum of two LinCombSpecies objects."""
        other = ensure_lincomb(other)
        return LinCombSpecies([*self.factor_species_list, *other.factor_species_list])

    def __sub__(self, other: SpeciesLike) -> LinCombSpecies:
        """Return LinCombSpecies object representing the difference of two LinCombSpecies objects."""
        return self + (-ensure_lincomb(other))

    def __eq__(self, value: Union[float, list[float], ndarray]) -> ConstraintLike:
        """Return Constraint object representing the equality of a LinCombSpecies object and a numerical value."""
        if isinstance(value, Union[list, ndarray]):
            return ConstraintSweep(self, value)
        return Constraint(self, value)

    def __rshift__(self, other: SpeciesLike) -> Reaction:
        """Return Reaction object representing a chemical reaction."""
        return Reaction(self, other)

    def print(self, return_str: bool = False) -> str:
        """Print the specifics of the LinCombSpecies object."""
        mystr = ""
        for fs in self.factor_species_list:
            if fs.factor > 0:
                mystr += "+"
            mystr += f"{fs.factor}*species_{fs.species_id} "
        if return_str:
            return mystr
        print(mystr)  # noqa: T201
        return ""


def ensure_lincomb(obj: SpeciesLike) -> LinCombSpecies:
    """Auxiliary function to ensure that a LinCombSpecies object is being used."""
    if isinstance(obj, Species):
        return LinCombSpecies([FactorSpecies(obj, 1.0)])
    if isinstance(obj, FactorSpecies):
        return LinCombSpecies([obj])
    return obj


# Type aliases
SpeciesLike = Union[Species, FactorSpecies, LinCombSpecies]
