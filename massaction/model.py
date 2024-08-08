"""Module for defining chemical species and reactions."""

from __future__ import annotations

import numpy as np
from scipy.optimize import root

from typing import Union


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
        return FactorSpecies(self.species, -1.0)

    def __mul__(self, factor: float) -> FactorSpecies:
        """Return FactorSpecies object representing the product of a species and a numerical factor."""
        return FactorSpecies(self, float(factor))

    def __rmul__(self, factor: float) -> FactorSpecies:
        """Return FactorSpecies object representing the product of a species and a numerical factor."""
        return FactorSpecies(self, float(factor))

    def __add__(self, other: ChemicalSet) -> LinCombSpecies:
        """Return LinCombSpecies object representing the sum of two Species objects."""
        if isinstance(other, Species):
            return LinCombSpecies([FactorSpecies(self, 1.0), FactorSpecies(other, 1.0)])
        if isinstance(other, FactorSpecies):
            return LinCombSpecies([FactorSpecies(self, 1.0), other])
        return LinCombSpecies([FactorSpecies(self, 1.0), *other.factor_species_list])

    def __sub__(self, other: Union[Species, FactorSpecies]) -> LinCombSpecies:
        """Return LinCombSpecies object representing the difference of two Species objects."""
        if isinstance(other, Species):
            return LinCombSpecies(
                [FactorSpecies(self, 1.0), FactorSpecies(other, -1.0)]
            )
        return LinCombSpecies([FactorSpecies(self, 1.0), -1.0 * other])

    def __eq__(self, value: float) -> Constraint:
        """Return Constraint object representing the equality of a Species object and a LinCombSpecies object."""
        return Constraint(self, value)

    def __rshift__(self, other: ChemicalSet) -> Reaction:
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

    def __add__(self, other: Union[Species, FactorSpecies]) -> LinCombSpecies:
        """Return LinCombSpecies object representing the sum of two FactorSpecies objects."""
        if isinstance(other, Species):
            return LinCombSpecies([self, FactorSpecies(other, 1.0)])
        return LinCombSpecies([self, other])

    def __radd__(self, other: Union[Species, FactorSpecies]) -> LinCombSpecies:
        """Return LinCombSpecies object representing the sum of a Species object and a FactorSpecies object."""
        return self.__add__(other)

    def __sub__(self, other: Union[Species, FactorSpecies]) -> LinCombSpecies:
        """Return LinCombSpecies object representing the difference of two FactorSpecies objects."""
        if isinstance(other, Species):
            return LinCombSpecies([self, FactorSpecies(other, -1.0)])
        return LinCombSpecies([self, -1.0 * other])

    def __eq__(self, value: float) -> Constraint:
        """Return Constraint object representing the equality of a FactorSpecies object and a numerical value."""
        return Constraint(self, value)

    def __rshift__(self, other: ChemicalSet) -> Reaction:
        """Return Reaction object representing a chemical reaction."""
        return Reaction(self, other)


class LinCombSpecies:
    """Class for modeling the linear combination of chemical species."""

    def __init__(self, factor_species_list: list) -> None:
        """Instantiate LinCombSpecies object.

        :param factor_species_list: List of FactorSpecies objects.
        """
        self.factor_species_list = factor_species_list

    def __add__(self, other: ChemicalSet) -> LinCombSpecies:
        """Return LinCombSpecies object representing the sum of two LinCombSpecies objects."""
        if isinstance(other, Species):
            return LinCombSpecies(
                [*self.factor_species_list, FactorSpecies(other, 1.0)]
            )
        if isinstance(other, FactorSpecies):
            return LinCombSpecies([*self.factor_species_list, other])
        return LinCombSpecies([*self.factor_species_list, *other.factor_species_list])

    def __radd__(self, other: ChemicalSet) -> LinCombSpecies:
        """Return LinCombSpecies object representing the sum of a Species object and a LinCombSpecies object."""
        return self.__add__(other)

    def __sub__(self, other: ChemicalSet) -> LinCombSpecies:
        """Return LinCombSpecies object representing the difference of two LinCombSpecies objects."""
        if isinstance(other, Species):
            return LinCombSpecies(
                [*self.factor_species_list, FactorSpecies(other, -1.0)]
            )
        if isinstance(other, FactorSpecies):
            return LinCombSpecies([*self.factor_species_list, -1.0 * other])
        return LinCombSpecies(
            [
                *self.factor_species_list,
                *[-1.0 * fs for fs in other.factor_species_list],
            ]
        )

    def __rsub__(self, other: Union[FactorSpecies, Species]) -> LinCombSpecies:
        """Return LinCombSpecies object representing the difference of a Species object and a LinCombSpecies object."""
        if isinstance(other, Species):
            return LinCombSpecies(
                [
                    *[-1.0 * fs for fs in self.factor_species_list],
                    FactorSpecies(other, 1.0),
                ]
            )
        if isinstance(other, FactorSpecies):
            return LinCombSpecies(
                [*[-1.0 * fs for fs in self.factor_species_list], other]
            )
        return LinCombSpecies(
            [
                *[-1.0 * fs for fs in other.factor_species_list],
                *other.factor_species_list,
            ]
        )

    def __eq__(self, value: float) -> Constraint:
        """Return Constraint object representing the equality of a LinCombSpecies object and a numerical value."""
        return Constraint(self, value)

    def __rshift__(self, other: ChemicalSet) -> Reaction:
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


def ensure_lincomb(obj: ChemicalSet) -> LinCombSpecies:
    """Auxiliary function to ensure that a LinCombSpecies object is being used."""
    if isinstance(obj, Species):
        return LinCombSpecies([FactorSpecies(obj, 1.0)])
    if isinstance(obj, FactorSpecies):
        return LinCombSpecies([obj])
    return obj


class Constraint:
    """Class for modeling constraints on chemical species."""

    def __init__(self, lincomb: ChemicalSet, value: float) -> None:
        """Instantiate Constraint object.

        :param lincomb: LinCombSpecies object.
        :param value: Value of the constraint.
        """
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


class Reaction:
    """Class for modeling chemical reactions.

    :param reactants: LinCombSpecies object representing the reactants.
    :param products: LinCombSpecies object representing the products.
    """

    def __init__(self, reactants: ChemicalSet, products: ChemicalSet) -> None:
        """Instantiate Reaction object."""
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


# Type aliases
ChemicalSet = Union[Species, FactorSpecies, LinCombSpecies]
