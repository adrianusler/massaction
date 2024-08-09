[![PyPI version](https://badge.fury.io/py/massaction.svg)](https://badge.fury.io/py/massaction)
[![CI Status](https://github.com/adrianusler/massaction/actions/workflows/test.yml/badge.svg)](https://github.com/adrianusler/massaction/actions/workflows/test.yml)
[![code coverage](https://img.shields.io/codecov/c/gh/adrianusler/massaction)](https://codecov.io/gh/adrianusler/massaction)

# massaction
A python package for the numerical solution of mass-action laws with constraints for the description of chemical reactions.

## Installation
The easiest way to install **massaction** is from [PyPi](https://pypi.org/project/massaction/) using pip.
```bash
pip install massaction
```

## Usage
The major perk of **massaction** is the simple syntax in which one can set up the description of a thermodynamic equilibrium. For instance, the mass-action law for the oxyhydrogen reaction may be set up in a couple of code lines:
```python
from massaction.model import ChemModel

mymodel = ChemModel(3) # set up a model with 3 chemical species
h2o, h2, o2 = mymodel.get_all_species() # give names to the species objects

# first set up the reaction equation, then define constraints to ensure atom balance
reaction = h2 + o2 >> 2*h2o
ln_equilibrium_constant = 20. # arbitrary units
constraint_hydrogen = 2*h2 + 2*h2o == 1.0 # arbitrary units
constraint_oxygen = h2o + 2*o2 == 10.0 # arbitrary units

# now solve the system of equations to obtain an array with the natural logarithm of the concentrations
ln_concentrations = model.solve( [reaction], [ln_equilibrium_constant], [constraint_hydrogen, constraint_oxygen] )
```
The resulting array `ln_concentrations` is `array([ -0.69314718, -22.94443898,   1.55814462])`, which contains the entries $\ln\left(c(\mathrm{H_2O})\right)$, $\ln\left(c(\mathrm{H_2})\right)$, and $\ln\left(c(\mathrm{O_2})\right)$, respectively.
