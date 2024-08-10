from massaction.model import ChemModel
from massaction.constraint import Constraint

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
        assert constraint.lincomb == lincomb
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
        assert constraint.lincomb == lincomb
        assert constraint.value == 10.7
        lnconc = np.array([3.14, 3.0, 2.9])
        conc = np.exp(lnconc)
        constraint_eval_pos = 3.4 * conc[0]
        constraint_eval_neg = 0.2 * conc[1] + conc[2] + 10.7
        assert np.isclose(
            constraint.eval(lnconc),
            np.log(constraint_eval_pos) - np.log(constraint_eval_neg),
        )
