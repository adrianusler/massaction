from massaction.model import ChemModel
import math  # standard library for mathematical functions

# DEFINE SIMPLE ACID DISSOCIATION MODEL
# HA <=> H+ + A-
# HA: undissociated acid, A-: acid anion, H+: proton
mymodel = ChemModel(3)
ha, a, h = mymodel.get_all_species()

# SET SIMULATION PARAMETERS
# - c0_ha: initial concentration of undissolved acid
# - pKa: acid dissociation constant
c0_ha = 0.1  # [mol/L]
pKa = 4.75
# In the current version, massaction requires the
# natural logarithms of the reactions' equilibrium constants
lnKa = math.log(10 ** (-pKa))

# DEFINE CONSTRAINTS
# 1.) Use the initial concentration of HA as a mathematical constraint
constraint1 = ha + a == c0_ha
# 2.) Electroneutrality constraint: the total charge is zero;
#     the total charge is the sum of the charges of all species.
#     A- and H+ are the only charged species in the system.
#     So the constraint is: [H+] - [A-] = 0
constraint2 = h - a == 0

# DEFINE REACTION EQUATIONS
# 1.) HA <=> H+ + A-
reaction1 = ha >> h + a

# SOLVE THE MODEL FOR EQUILIBRIUM CONCENTRATIONS
# In the current version, massaction returns the
# natural logarithms of the equilibrium concentrations.
lnc_equilibrium = mymodel.solve(
    [reaction1],
    [lnKa],
    [constraint1, constraint2],
)
c_equilibrium = [math.exp(lnc) for lnc in lnc_equilibrium]

# PRINT RESULTS
print(
    f"Initial concentration of HA: {c0_ha:.2e} mol/L\n"
    "Equilibrium concentrations:\n"
    f"  HA: {c_equilibrium[0]:.2e} mol/L\n"
    f"  A-: {c_equilibrium[1]:.2e} mol/L\n"
    f"  H+: {c_equilibrium[2]:.2e} mol/L"
)
