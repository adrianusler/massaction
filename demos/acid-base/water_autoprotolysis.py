from massaction.model import ChemModel, nil
import math  # standard library for mathematical functions

# DEFINE SIMPLE AUTOPROTOLYSIS MODEL
# nil <=> H+ + OH-
# nil: "nothing"; actually, two water molecules react.
# The concentration of water is so high that it is considered constant.
# The concentration of water is, therefore, incorporated into the equilibrium constant,
# i.e., the autoprotolysis constant, Kw = [H+][OH-]
mymodel = ChemModel(2)
h, oh = mymodel.get_all_species()

# SET SIMULATION PARAMETERS
# - c0_ha: initial concentration of undissolved acid
# - Kw: ionic product of water
Kw = 1e-14  # [mol^2/L^2], value at 25°C
# In the current version, massaction requires the
# natural logarithms of the reactions' equilibrium constants
lnKw = math.log(Kw)

# DEFINE CONSTRAINTS
#     Electroneutrality constraint: the total charge is zero;
#     the total charge is the sum of the charges of all species.
#     OH- and H+ are the only charged species in the system.
#     So the constraint is: [H+] - [OH-] = 0
constraint = h - oh == 0

# DEFINE REACTION EQUATIONS
reaction = nil >> h + oh

# SOLVE THE MODEL FOR EQUILIBRIUM CONCENTRATIONS
# In the current version, massaction returns the
# natural logarithms of the equilibrium concentrations.
lnc_equilibrium = mymodel.solve(
    [reaction],
    [lnKw],
    [constraint],
)
c_equilibrium = [math.exp(lnc) for lnc in lnc_equilibrium]

# PRINT RESULTS
print(
    "Equilibrium concentrations:\n"
    f"   H+: {c_equilibrium[0]:.2e} mol/L\n"
    f"  OH-: {c_equilibrium[1]:.2e} mol/L"
)
