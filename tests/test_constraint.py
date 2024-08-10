from massaction.model import ChemModel
from massaction.constraint import Constraint, ConstraintSweep
from massaction.species import ensure_lincomb

from unittest import TestCase
import numpy as np


class TestConstraint(TestCase):
    """Test the Constraint class."""

    def setUp(self):
        self.model3 = ChemModel(3)

    def test_constraint(self):
        h2o, h2, o2 = self.model3.get_all_species()

        lincomb = 3.4 * h2o - 0.2 * h2 - o2
        constraint = lincomb == -1.3
        assert isinstance(constraint, Constraint)
        # assert constraint.lincomb == lincomb
        assert constraint.value == -1.3
        lnconc = np.array([23.0, -1.0, 4.3])
        conc = np.exp(lnconc)
        constraint_eval_pos = 3.4 * conc[0] + 1.3
        constraint_eval_neg = 0.2 * conc[1] + conc[2]
        assert np.isclose(
            constraint.eval(lnconc),
            np.log(constraint_eval_pos) - np.log(constraint_eval_neg),
        )

        constraint = lincomb == 10.7
        assert isinstance(constraint, Constraint)
        assert constraint.lincomb is lincomb
        assert constraint.value == 10.7
        lnconc = np.array([3.14, 3.0, 2.9])
        conc = np.exp(lnconc)
        constraint_eval_pos = 3.4 * conc[0]
        constraint_eval_neg = 0.2 * conc[1] + conc[2] + 10.7
        assert np.isclose(
            constraint.eval(lnconc),
            np.log(constraint_eval_pos) - np.log(constraint_eval_neg),
        )


class TestConstraintSweep(TestCase):
    """Test the ConstraintSweep class."""

    def setUp(self):
        self.model3 = ChemModel(3)
        self.h2o, self.h2, self.o2 = self.model3.get_all_species()

    def test_constraint_sweep(self):
        def cs_assertions(cs, lincomb, values):
            assert isinstance(cs, ConstraintSweep)
            for fs1, fs2 in zip(
                cs.lincomb.factor_species_list, lincomb.factor_species_list
            ):
                assert fs1.factor == fs2.factor
                assert fs1.species is fs2.species
            assert cs.num_values == len(values)
            assert cs.current_id == 0
            assert cs.value == values[0]
            for i in range(len(values)):
                cs.set_current_id(i)
                assert cs.current_id == i
                assert cs.values[i] == values[i]
                assert cs.value == values[i]
            # test cs.eval ?

        h2o, h2, o2 = self.h2o, self.h2, self.o2
        # test Species.__eq__
        cs_assertions(h2o == [1.0, 2.0, 3.0], ensure_lincomb(h2o), [1.0, 2.0, 3.0])
        cs_assertions(
            o2 == np.linspace(0.4, 0.7, 5), ensure_lincomb(o2), np.linspace(0.4, 0.7, 5)
        )
        # test FactorSpecies.__eq__
        cs_assertions(
            7 * h2o == [6.0, 1.0, 8.0], ensure_lincomb(7 * h2o), [6.0, 1.0, 8.0]
        )
        cs_assertions(
            2 * o2 == np.linspace(-2.0, 9, 3),
            ensure_lincomb(2 * o2),
            np.linspace(-2.0, 9, 3),
        )
        # test LinCombSpecies.__eq__
        cs_assertions(h2o - h2 + o2 == [2.0, 5.0, 4.0], h2o - h2 + o2, [2.0, 5.0, 4.0])
        cs_assertions(
            3 * o2 - h2 == np.linspace(12.0, -11.0, 4),
            3 * o2 - h2,
            np.linspace(12, -11.0, 4),
        )
