# Computing consumption fractions from competition between species with optimal parameters from 'main_single_species.py'
# data and competitor
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import time

from model_equations import a_b, ap_bp, lag_min, delta_max
from simulation_functions import run_competition_in_parallel

np.random.seed(18)

font = {'family': 'Times New Roman',
        'weight': 'normal',
        'size': 20}
mpl.rc('font', **font)

#lag_cmap = mpl.cm.get_cmap('viridis')
#del_cmap = mpl.cm.get_cmap('plasma')
constant_index = {'T0':0, 'Tab':1}



###########################
## Simulation parameters ##
###########################
bac_res = 12                           # resolution in bacterial parameters
ab_res = 11                            # resolution in antibiotic parameters
t_res = 10                                  # resolution in time array
tot_cycles = 20
repetitions = 1                             # number of repetitions for ensemble average

data = 'new'                                # 'new' - generate and plot new data. 'old' plot old data
save_fig = False
save_data = False



###########################
## Antibiotic parameters ##
###########################
# for chosing which time parameter to keep constant.
ic = constant_index['Tab']                      # 'T0' or 'Tab'
T_const = 8                                     # value of the constant parameter

# defining parameter arrays
T_max = [12, 24]                                # upper bounds on meningful values for T0 and Tab
p_arr = np.linspace(0, 1, ab_res)               # probability of antibiotics
T_arr = np.linspace(0, T_max[1-ic], ab_res)     # time array

T_labels = ['T0', 'Tab']
T_values = [T_const, T_arr]
T0 = T_values[ic]
Tab = T_values[1-ic]



##########################
## Bacterial parameters ##
##########################
# empty arrays parameters
lag = np.ones([2, bac_res, bac_res]) * lag_min
delta = np.zeros([2, bac_res, bac_res])

# competitor parameters
I = np.ones(bac_res)
lag[1] += np.outer(I, np.linspace(0, max(T0 + Tab), bac_res))
delta[1] += np.outer(np.linspace(0, delta_max, bac_res), I)

# transforming parameters to a-b space
a, b = a_b(lag, delta)
ap, bp = ap_bp(lag, delta)



########################
## Running Simulation ##
########################
if data == 'new':
    tic = time.time()
    bac_args = [lag, delta, a, b, ap, bp]
    ab_args = [p_arr, T0, Tab]
    sim_args = [ab_res, bac_res, t_res, tot_cycles, repetitions, save_data]

    if __name__ == '__main__':
        S_frac, lag_opt, del_opt = run_competition_in_parallel(bac_args, ab_args, sim_args)

    # saving
    if save_data:
        np.savetxt("../../../data/model2/competition_Sfrac-"+T_labels[ic]+str(T_const), S_frac)
        np.savetxt("../../../data/model2/competition_lag-"+T_labels[ic]+str(T_const), lag_opt)
        np.savetxt("../../../data/model2/competition_delta-"+T_labels[ic]+str(T_const), del_opt)

    toc = time.time()
    print(f'\nTime: {(toc - tic) / 3600}h\n')


###################
## Plotting data ##
###################
if data == 'old':
    S_frac = np.loadtxt('../../../data/model2/competition_Sfrac-'+T_labels[ic]+str(T_const))[:,:-1]
    lag_opt = np.loadtxt('../../../data/model2/competition_lag-'+T_labels[ic]+str(T_const))[:,:-1]
    del_opt = np.loadtxt('../../../data/model2/competition_delta-'+T_labels[ic]+str(T_const))[:,:-1]


T = np.outer(np.ones(ab_res), T0) + np.outer(np.ones(ab_res), Tab)
T[0, 0] = 1
time_labels = [r'$T_0$', r'$T_{AB}$']
title = time_labels[ic] + ' = ' + str(T_const)

fig, ax = plt.subplots(1, 2, figsize=(15, 5.5), sharey=True)
# fig.suptitle(title)

# plotting lag time
ax[0].set(xlabel=time_labels[1-ic] + r'$~[h]$', ylabel=r"$p$")#, title=r"$\lambda^*$")
im0 = ax[0].imshow(lag_opt / T, origin="lower", cmap=lag_cmap, aspect="auto", vmin=0, extent=[0, T_arr.max(), 0, 1])
cbar = fig.colorbar(im0, ax=ax[0], aspect=20, anchor=(-.1, 0.5))
cbar.set_label(r"$\lambda^* / ~T$", rotation=0, labelpad=30)

# plotting delta
ax[1].set(xlabel=time_labels[1-ic] + r'$~[h]$', ylabel=r"$p$")#, title=r"$\delta^*$")
im1 = ax[1].imshow(del_opt, origin="lower", cmap=del_cmap, aspect="auto", vmin=0, vmax=0.075, extent=[0, T_arr.max(), 0, 1])
cbar = fig.colorbar(im1, ax=ax[1], aspect=20, anchor=(-.1, 0.5))
cbar.set_label(r"$\delta^* [h^{-1}]$", rotation=0, labelpad=50)

# ax[2].set(xlabel=time_labels[1-ic] + r'$~[h]$', ylabel=r"$p$")#, title="Consumption fraction")
# im2 = ax[2].imshow(S_frac, origin="lower", aspect="auto", extent=[0, max(T_arr), 0, 1])
# cbar = fig.colorbar(im2, ax=ax[2], aspect=20, anchor=(-.1, 0.5))
# cbar.set_label("Consumption fraction", rotation=0, labelpad=50)

fig.tight_layout()
fig.show()

if save_fig:
    fig.savefig("../../../figures/model2/competitor_parameters-"+T_labels[ic]+str(T_const)+".png", dpi=100)

3