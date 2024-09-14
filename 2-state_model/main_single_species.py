# computing optimal parameter of single species with restricted nutrients
import sys
import numpy as np
import time

sys.path.append("src")
from model_equations import a_b, ap_bp, lag_min, delta_max
from simulation_functions import compute_optimal_parameters, analytical_fitness

constant_index = {'T0':0, 'Tab':1}

###########################
## Simulation parameters ##
###########################
save_data = False
data = 'old'                                    # 'new' - generate and plot new data. 'old' plot old data
resolution = 'high'                              # resolution of antibiotic parameters. 'low' or 'high'

if resolution == 'low': ab_res = 26 #101               # for generating data used in competition simulations
elif resolution == 'high': ab_res = 400            # for creating high resolution plots
bac_res = 400                                   # resolution in bacterial parameters



###########################
## Antibiotic parameters ##
###########################
# for chosing which time parameter to keep constant.
ic = constant_index['Tab']                      # 'T0' or 'Tab'
T_const = 10                                    # value of the constant parameter

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
bac_args = [lag, delta, a, b, ap, bp]
ab_args = [p_arr, T0, Tab]
tic = time.time()
lag_opt, del_opt, fitness = compute_optimal_parameters(bac_args, ab_args)
toc = time.time()
print(f'Time: {toc - tic}s')

# saving
np.savetxt(f"data/half_dilution_ratio/optimal_lag-{T_labels[ic]}{T_const}.txt", lag_opt)
np.savetxt(f"data/half_dilution_ratio/optimal_delta-{T_labels[ic]}{T_const}.txt", del_opt)
np.savetxt(f"data/half_dilution_ratio/optimal_fitness-{T_labels[ic]}{T_const}.txt", fitness)

