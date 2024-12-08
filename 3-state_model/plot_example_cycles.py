import sys
import argparse
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

sys.path.append("src")
from differential_equations import S0, f, λ_min
from differential_equations import ode_grow, ode_kill

mpl.rcParams["font.size"] = "12"

parser = argparse.ArgumentParser(description='Competition between N species for tot_cycles cycles.')
parser.add_argument('T0',  type=float, help='application time of antibiotics')
parser.add_argument('Tab', type=float, help='duration of antibiotics')
args = parser.parse_args()

T0  = args.T0
Tab = args.Tab
T   = T0 + Tab
T_cycle = 40

# index of substrate
iS = 6

p = 0.6
r_arr = [0, 1, 1, 0, 1]
bac_args=[np.array([T, λ_min]), np.array([λ_min, 3]), np.array([0, 0.04])]


popu_0 = [f*S0, f*S0, 0, 0, 0, 0, S0]
popu_t = [[], [], [], [], [], [], []]
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
    popu_0 = f * np.array([popu_t[0][ic][-1] + popu_t[2][ic][-1] + popu_t[4][ic][-1], popu_t[1][ic][-1] + popu_t[3][ic][-1] + popu_t[5][ic][-1], 0, 0, 0, 0, S0 / f])
    


##### PREPARE ARRAYS ####
# concatenating arrays for plotting
for i in range(iS):
    popu_t[i] = np.concatenate(popu_t[i])

species1 = popu_t[0] + popu_t[2] + popu_t[4]
species2 = popu_t[1] + popu_t[3] + popu_t[5]

S = np.concatenate(popu_t[6])
time = np.concatenate(tot_time)



##### PLOTTING #####
sns.set_theme(style='ticks', palette='deep', font_scale=1.2)
fig, ax = plt.subplots(1, 1, figsize=(6.75, 2.7))
ax.fill_between(time, S, 0, color=sns.color_palette('deep')[2], alpha=0.3, label="Substrate")

ax.plot(time, species1, lw=2, color=sns.color_palette('deep')[0], label="Species 1")
ax.plot(time, species2, lw=2, color=sns.color_palette('deep')[3], label="Species 2")
ax.set_xticks(np.arange(0, max(time), max(time)/5) + max(time)/10)
ax.set_xticklabels(['1', '2', '3', '4', '5'])

ax.set(xlabel="Cycle", ylabel="Population")
ax.set(yscale='log',  ylim=[10**(0), 10**10])
sns.despine()

fig.tight_layout(rect=[0, 0, 1, 0.85])
ax.legend(loc='upper center',
          bbox_to_anchor=(0.5, 1.45),
          ncol=3, 
          frameon=False)
fig.savefig("figs/example_3state.png")
