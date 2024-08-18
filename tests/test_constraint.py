from massaction.model import ChemModel
from massaction.constraint import Constraint, ConstraintSweep
from massaction.species import ensure_lincomb

from unittest import TestCase
import numpy as np

import pytest


class TestConstraint(TestCase):
    """Test the Constraint class."""

    def setUp(self):
        self.model3 = ChemModel(3)

    def test_constraint(self):
        def cstr_assertions(cstr, lincomb, value):
            lnc = np.array([4.7, -6.2, 7.3])
            c = np.exp(lnc)
            cstr_eval_pos = 0.0
            cstr_eval_neg = 0.0
            assert isinstance(cstr, Constraint)
            for fs1, fs2 in zip(
                cstr.lincomb.factor_species_list, lincomb.factor_species_list
            ):
                assert fs1.factor == fs2.factor
                assert fs1.species is fs2.species
                if fs1.factor > 0:
                    cstr_eval_pos += fs1.factor * c[fs1.species.species_id]
                else:
                    cstr_eval_neg -= fs1.factor * c[fs1.species.species_id]
            if value > 0:
                cstr_eval_neg += value
            else:
                cstr_eval_pos -= value
            assert np.isclose(
                cstr.eval(lnc),
                np.log(cstr_eval_pos) - np.log(cstr_eval_neg),
            )
            assert cstr.value == value
            assert cstr.num_values == 1  # by definition, since it is not a sweep

        h2o, h2, o2 = self.model3.get_all_species()

        # test LinCombSpecies.__eq__
        cstr = (3.4 * h2o - 0.2 * h2 - o2) == -1.3
        cstr_assertions(
            cstr,
            ensure_lincomb(3.4 * h2o - 0.2 * h2 - o2),
            -1.3,
        )
        assert cstr.check_reservoir() == (False, -1)
        assert cstr.reservoir_value == 0.0
        # test FactorSpecies.__eq__
        cstr_reservoir = 2 * h2o == 10.7
        cstr_assertions(2 * h2o == 10.7, ensure_lincomb(2 * h2o), 10.7)
        assert cstr_reservoir.check_reservoir() == (True, 0)
        assert cstr_reservoir.reservoir_value == 5.35
        # test Species.__eq__
        cstr_reservoir = h2 == 4.5
        cstr_assertions(cstr_reservoir, ensure_lincomb(h2), 4.5)
        assert cstr_reservoir.check_reservoir() == (True, 1)
        assert cstr_reservoir.reservoir_value == 4.5

    @pytest.fixture(autouse=True)
    def capsys(self, capsys):
        self.capsys = capsys

    def test_print(self):
        h2o, h2, o2 = self.model3.get_all_species()
        # test Species.__eq__
        (h2 == 4.5).print()
        out, _ = self.capsys.readouterr()
        assert out == "+1.0*species_1  == 4.5\n"
        # test FactorSpecies.__eq__
        (2 * h2o == 10.7).print()
        out, _ = self.capsys.readouterr()
        assert out == "+2.0*species_0  == 10.7\n"
        # test LinCombSpecies.__eq__
        (h2 - h2o == -1.3).print()
        out, _ = self.capsys.readouterr()
        assert out == "+1.0*species_1 -1.0*species_0  == -1.3\n"


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
