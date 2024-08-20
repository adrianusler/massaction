from massaction.model import ChemModel, nil
import math  # standard library for mathematical functions

# AH <=> H+ + A-
# AH: undissociated acid; A-: acid anion; H+: proton
# nil <=> H+ + OH-
# nil: "nothing". Actually, two water molecules react.
# The concentration of water is so high that it is considered constant.
# The concentration of water is, therefore, incorporated into the equilibrium constant,
# i.e., the autoprotolysis constant, Kw = [H+][OH-]
mymodel = ChemModel(4)
h, oh, ha, a = mymodel.get_all_species()

# SET SIMULATION PARAMETERS
# - c0_ha: initial concentration of undissolved acid
# - Kw: ionic product of water
Kw = 1e-14  # [mol^2/L^2], value at 25°C
# In the current version, massaction requires the
# natural logarithms of the reactions' equilibrium constants
lnKw = math.log(Kw)
# - c0_ha: initial concentration of undissolved acid
# - pKa: acid dissociation constant
c0_ha = 0.1  # [mol/L]
pKa = 4.75
# In the current version, massaction requires the
# natural logarithms of the reactions' equilibrium constants
lnKa = math.log(10 ** (-pKa))


# DEFINE CONSTRAINTS
#     Electroneutrality constraint: the total charge is zero;
#     the total charge is the sum of the charges of all species.
#     OH- and H+ are the only charged species in the system.
#     So the constraint is: [H+] - [OH-] = 0
constraint_elec = h - a - oh == 0
#     Constraint for the acid rest "A"
constraint_acid = ha + a == c0_ha

# DEFINE REACTION EQUATIONS
# 1.) HA <=> H+ + A-
reaction_acid = ha >> a + h
# 2.) nil <=> H+ + OH-
reaction_water = nil >> h + oh

# SOLVE THE MODEL FOR EQUILIBRIUM CONCENTRATIONS
# In the current version, massaction returns the
# natural logarithms of the equilibrium concentrations.
lnc_equilibrium = mymodel.solve(
    [reaction_acid, reaction_water],
    [lnKa, lnKw],
    [constraint_elec, constraint_acid],
)

c_h_eq, c_oh_eq, c_ha_eq, c_a_eq = [math.exp(lnc) for lnc in lnc_equilibrium]

# PRINT RESULTS
print(
    "Equilibrium concentrations:\n"
    f"   H+: {c_h_eq:.2e} mol/L\n"
    f"  OH-: {c_oh_eq:.2e} mol/L\n"
    f"   HA: {c_a_eq:.2e} mol/L\n"
    f"   A-: {c_ha_eq:.2e} mol/L\n"
)
