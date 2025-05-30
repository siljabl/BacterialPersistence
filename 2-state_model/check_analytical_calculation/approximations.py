import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../src")
from analytical_calculations import compute_a_and_b, compute_ap_and_bp, solve_constants
from differential_equations import lag_min, delta_max

res = 100
T_max = 6
lag_arr   = np.linspace(0, T_max, res) + lag_min
delta_arr = np.linspace(0, delta_max, res)
T0 = 0
T = T_max - T0

lag   = np.outer(np.ones(res), lag_arr)
delta = np.outer(delta_arr, np.ones(res))

a,  b  = compute_a_and_b(lag, delta)
ap, bp = compute_ap_and_bp(lag, delta)

eq_params = {'a':a, 'b':b, 'ap':ap, 'bp':bp}
ab_params = {'T0':T0, 'T':T}

#C1, C2 = solve_constants(eq_params, ab_params, stage="pre")
D1, D2 = solve_constants(eq_params, ab_params, stage="ab")

c = 1/100

plt.imshow((b+a)*(bp-c)*(ap-c)/((a-c)*(b+ap)*(b+bp)), vmin=-10)
plt.colorbar()
plt.xlabel(r"$\lambda$")
plt.ylabel(r"$\delta$")
plt.savefig("a_minus.png")

# T = np.linspace(lag_min, T_max, 100)
# value = np.sqrt(1 + 4*T) / T

# plt.plot(lag_arr, a-1/lag_arr)
# plt.plot(lag_arr, a-ap)
#plt.plot(lag_arr, a-bp)

# a,  b  = compute_a_and_b(lag_arr, 0.001)
# plt.plot(lag_arr, 1/a)
# plt.savefig("a_minus.png")

print(delta_max)