from massaction.model import ChemModel
import numpy as np

mymodel = ChemModel(3)
nil, vm, vo = mymodel.get_all_species()

# Define constraints in terms of molar fractions (particles per formula unit of crystal)
constraint_nil = nil == 1.0  # This is a dummy variable
constraint_electroneutrality = 2 * vo - 4 * vm == 0.0

# Define the reaction equation
reaction = nil >> vm + 2 * vo
lnK = -25.0

lnc = mymodel.solve([reaction], [lnK], [constraint_nil, constraint_electroneutrality])
c = np.exp(lnc)
print(
    "Results:"
    f"\n\t[nil] = {c[0]:.2e}"
    f"\n\t[vm] = {c[1]:.2e} (i.e., {c[1]*100:1.3f}% of M sites)"
    f"\n\t[vo] = {c[2]:.2e} (i.e., {.5*c[2]*100:1.3f}% of O sites)"
)
