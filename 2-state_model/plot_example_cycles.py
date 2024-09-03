import sys
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

sys.path.append("src")
from model_equations import S0, f, lag_min
from model_equations import ode_grow, ode_kill

parser = argparse.ArgumentParser(description='Competition between N species for tot_cycles cycles.')
parser.add_argument('T0',  type=float, help='application time of antibiotics')
parser.add_argument('Tab', type=float, help='duration of antibiotics')
args = parser.parse_args()

T0  = args.T0
Tab = args.Tab
T   = T0 + Tab
T_cycle = 40

# index of substrate
iS = 4

p = 0.6
r_arr = [0, 1, 1, 0, 1]
bac_args=[np.array([T, 4]), np.array([lag_min, 0.01])]


popu_0 = [f*S0, f*S0, 0, 0, S0]
popu_t = [[], [], [], [], []]
cycle_time, tot_time = [], []

for ic in range(5):
    r = r_arr[ic]

    if r < p:
        # solving ODEs
        sol_cycle1 = solve_ivp(ode_grow, [0, T0],      popu_0,              args=bac_args, max_step=0.1)
        sol_cycle2 = solve_ivp(ode_kill, [T0, T],      sol_cycle1.y[:, -1], args=bac_args, max_step=0.1)
        sol_cycle3 = solve_ivp(ode_grow, [T, T_cycle], sol_cycle2.y[:, -1], args=bac_args, max_step=0.1)

        # assembling arrays
        t_tmp    = np.concatenate([sol_cycle1.t,   sol_cycle2.t,   sol_cycle3.t])
        popu_tmp = np.concatenate([sol_cycle1.y.T, sol_cycle2.y.T, sol_cycle3.y.T]).T

        t = max(sol_cycle3.t[sol_cycle3.y[iS] > 0])
        # time.append(t + ic * T_cycle)

    else:
        # solving ODE
        sol_cycle = solve_ivp(ode_grow, [0, T_cycle], popu_0, args=bac_args, max_step=0.1)
        t_tmp    = sol_cycle.t
        popu_tmp = sol_cycle.y

        t = max(sol_cycle.t[sol_cycle.y[iS] > 0])
        #cycle_time.append(t + ic * T_cycle)

    # saving data
    tot_time.append(t_tmp + ic * T_cycle)
    for i in range(iS+1):
        popu_t[i].append(popu_tmp[i])

    # preparing new cycle
    popu_0 = f * np.array([popu_t[0][ic][-1] + popu_t[2][ic][-1], popu_t[1][ic][-1] + popu_t[3][ic][-1], 0, 0, S0 / f])
    


##### PREPARE ARRAYS ####
# concatenating arrays for plotting
for i in range(iS):
    popu_t[i] = np.concatenate(popu_t[i])

species1 = popu_t[0] + popu_t[2]
species2 = popu_t[1] + popu_t[3]

S = np.concatenate(popu_t[iS])
time = np.concatenate(tot_time)



##### PLOTTING #####
#mpl.rcParams["text.usetex"] = True
mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.size"] = "12"

fig, ax = plt.subplots(1, 1, figsize=(6.5, 2.8))
ax.fill_between(time, S, 0, color="lightseagreen", alpha=0.3, label="substrate")

ax.plot(time, species1, lw=2, color="firebrick", label=rf"$(\lambda, \delta) = ({bac_args[0][0]:2.0f}, {bac_args[1][0]:2.2f})$")
ax.plot(time, species2, lw=2, color="royalblue", label=rf"$(\lambda, \delta) = ({bac_args[0][1]:2.0f}, {bac_args[1][1]:2.2f})$")

ax.set(xlabel="Time", ylabel="log(Population)")
ax.set(yscale='log',  ylim=[10**(-2), 10**10])

fig.tight_layout(rect=[0, 0.13, 1, 1])
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=3, fancybox=True, shadow=False)
fig.savefig("figs/example_2state.png")

# Chose illustrative bacterial parameters
# match size with text