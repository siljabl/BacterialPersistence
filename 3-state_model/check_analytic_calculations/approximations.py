import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../src")
from analytical_calculations import compute_a_and_b, compute_ap_and_bp, solve_constants
from differential_equations import λ_min, δ_max

res = 100
T_max = 36

T0 = 12
p_arr   = np.linspace(0, 1, res)
Tab_arr = np.linspace(0, T_max-12, res)


##########################
## Bacterial parameters ##
##########################
λ_arr = np.linspace(0, T_max, res) + λ_min
δ_arr = np.linspace(0, δ_max, res)

λd = np.outer(np.ones(res), np.outer(np.ones(res), λ_arr)).reshape(res, res, res) - λ_min / 1000    # avoid overflow by distinguishing λd and λr
λr = np.outer(np.outer(np.ones(res), λ_arr), np.ones(res)).reshape(res, res, res)
δ  = np.outer(np.outer(δ_arr, np.ones(res)), np.ones(res)).reshape(res, res, res)


###########################
## Simulation Parameters ##
###########################
a, b   = compute_a_and_b(λr, δ)
ap, bp = compute_ap_and_bp(λr, δ)
c = 1/λd

bac_params = {'λd':λd, 'λr':λr, 'δ':δ}
ab_params  = {'T0':T0, 'Tab':Tab_arr, 'T':T0+Tab_arr}
eq_params  = {'a':a, 'b':b, 'c':c, 'ap':ap, 'bp':bp}

B, A, C = solve_constants(eq_params, ab_params, stage="post")

t = 4
#plt.imshow(np.mean(abs(a-c), axis=0), vmax=1, vmin=0)
plt.imshow(abs(a-c)[:,0], vmax=1, vmin=0)
plt.colorbar()
plt.xlabel(r"$\lambda$")
plt.ylabel(r"$\delta$")
plt.savefig("a_minus.png")