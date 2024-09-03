# computing optimal parameter of single species with restricted nutrients
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import time

from model_equations import a_b, ap_bp, lag_min, delta_max
from simulation_functions import compute_optimal_parameters, analytical_fitness

font = {'family': 'Times New Roman',
        'weight': 'normal',
        'size': 20}
mpl.rc('font', **font)

lag_cmap = mpl.cm.get_cmap('viridis')
del_cmap = mpl.cm.get_cmap('plasma')
constant_index = {'T0':0, 'Tab':1}

###########################
## Simulation parameters ##
###########################
save_fig = False
save_data = False
data = 'old'                                    # 'new' - generate and plot new data. 'old' plot old data
resolution = 'high'                              # resolution of antibiotic parameters. 'low' or 'high'

if resolution == 'low': ab_res = 26 #101               # for generating data used in competition simulations
elif resolution == 'high': ab_res = 400            # for creating high resolution plots
bac_res = 25 #400                                   # resolution in bacterial parameters



###########################
## Antibiotic parameters ##
###########################
# for chosing which time parameter to keep constant.
ic = constant_index['T0']                      # 'T0' or 'Tab'
T_const = 5                                    # value of the constant parameter

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
lag_arr = np.linspace(0, T_const + T_arr.max() + 1, bac_res) + lag_min
delta_arr = np.linspace(0, delta_max, bac_res)

lag = np.outer(np.ones(bac_res), lag_arr)
delta = np.outer(delta_arr, np.ones(bac_res)) #+ 10**(-6)

a, b = a_b(lag, delta)
ap, bp = ap_bp(lag, delta)


########################
## Running Simulation ##
########################
if data == 'new':
    bac_args = [lag, delta, a, b, ap, bp]
    ab_args = [p_arr, T0, Tab]
    tic = time.time()
    lag_opt, del_opt, fitness = compute_optimal_parameters(bac_args, ab_args)
    toc = time.time()
    print(f'Time: {toc - tic}s')

    # saving
    if save_data:
        np.savetxt("../../../data/model2/"+resolution+"_resolution/optimal_lag-"+T_labels[ic]+str(T_const), lag_opt)
        np.savetxt("../../../data/model2/"+resolution+"_resolution/optimal_delta-"+T_labels[ic]+str(T_const), del_opt)
        np.savetxt("../../../data/model2/"+resolution+"_resolution/optimal_fitness-"+T_labels[ic]+str(T_const), fitness)


###################
## Plotting data ##
###################
# importing data 
if data == 'old':
    lag_opt = np.loadtxt('../../../data/model2/'+resolution+'_resolution/optimal_lag-' + T_labels[ic]+str(T_const))
    del_opt = np.loadtxt('../../../data/model2/'+resolution+'_resolution/optimal_delta-' + T_labels[ic]+str(T_const))
    fitness = np.loadtxt('../../../data/model2/'+resolution+'_resolution/optimal_fitness-' + T_labels[ic]+str(T_const))

# arrays for scaling imshow
T = np.outer(np.ones_like(Tab), Tab) + np.outer(np.ones_like(T0), T0)
p_mat = np.outer(p_arr, np.ones_like(p_arr))
T[0, 0] = 1
p_mat[0] = 1

# preparing strings
time_labels = [r'$T_0$', r'$T_{AB}$']
title = time_labels[ic] + ' = ' + str(T_const)

# setting up figure
fig, ax = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)
# fig.suptitle(title)

# plotting lag
ax[0].set(xlabel=time_labels[1-ic] + r'$~[h]$', ylabel=r"$p$")#, title=r"$\lambda^*$")
im0 = ax[0].imshow(lag_opt / T, origin="lower", cmap=lag_cmap, aspect="auto", vmin=0, extent=[T_arr.min(), T_arr.max(), 0, 1])
cbar = fig.colorbar(im0, ax=ax[0], aspect=20, anchor=(-.1, 0.5))
cbar.set_label(r"$\lambda^* / ~T$", rotation=0, labelpad=30)

# plotting delta
ax[1].set(xlabel=time_labels[1-ic] + r'$~[h]$', ylabel=r"$p$")#, title=r"$\delta^*$")
im1 = ax[1].imshow(del_opt, origin="lower", cmap=del_cmap, aspect="auto", vmin=0, vmax=0.075, extent=[T_arr.min(), T_arr.max(), 0, 1])
cbar = fig.colorbar(im1, ax=ax[1], aspect=20, anchor=(-.1, 0.5))
cbar.set_label(r"$\delta^* [h^{-1}]$", rotation=0, labelpad=50)

# plotting Ts
# ax[2].set(xlabel=time_labels[1-ic] + r'$~[h]$', ylabel=r"$p$")#, title=r"$T_S$")
# im2 = ax[2].imshow(fitness, origin="lower", aspect="auto", extent=[T_arr.min(), T_arr.max(), 0, 1])
# cbar = fig.colorbar(im2, ax=ax[2], aspect=20, anchor=(-.1, 0.5))
# cbar.set_label(r"$\langle T_S \rangle_p^{-1} ~[h^{-1}]$", rotation=0, labelpad=50)

# saving
fig.tight_layout()
plt.show()

if save_fig:
    fig.savefig("../../../figures/model2/optimal_parameters-"+T_labels[ic]+str(T_const)+".png", dpi=100)

