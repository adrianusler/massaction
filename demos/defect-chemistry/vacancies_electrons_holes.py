"""Demo for modelling the defect chemistry of an acceptor-doped oxide.
It involves the incorporation of the acceptor dopant and the formation of electron/hole pairs."""

from massaction.model import ChemModel, nil
import numpy as np
import matplotlib.pyplot as plt

# Define the model and then the species:
# - o2: oxygen gas
# - vo: oxygen vacancies
# - ox: oxygen ions
# - em: electrons
# - hp: holes
mymodel = ChemModel(5)
o2, vo, ox, em, hp = mymodel.get_all_species()

# Parameters
xdop = 0.05# acceptor-doping level

# Define the constraints
cstr_o2 = ( o2 == np.logspace(-35, 0, 100) )# variation of oxygen partial pressure
cstr_en = 2*vo + hp - em == xdop# electroneutrality constraint
cstr_ox = ox + vo == 2.

# Define the reactions
rct_incorp = .5*o2 + vo >> ox + 2*hp
rct_pair = nil >> em + hp

# Set values of the equilibrium constants
lnK_incorp = 5.
lnK_pair = -20.

# Solve the system
lnc = mymodel.solve(
    [rct_incorp, rct_pair],
    [lnK_incorp, lnK_pair],
    [cstr_o2, cstr_en, cstr_ox],
)
c_o2, c_vo, c_ox, c_em, c_hp = np.exp(lnc).T

# Plot the results
plt.plot(c_o2, c_em, color="C0", label=r'$\mathrm{e^\prime}$')
plt.plot(c_o2, c_hp, color="C1", label=r'$\mathrm{h^\bullet}$')
plt.plot(c_o2, c_vo, color="C3", label=r'$\mathrm{v_O^{\bullet\bullet}}$')
plt.xscale("log")
plt.yscale("log")
plt.xlabel("oxygen partial pressure (bar)")
plt.ylabel("defects per formula unit")
plt.legend()
plt.show()