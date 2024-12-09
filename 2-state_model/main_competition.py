# Computing consumption fractions from competition between species with optimal parameters from 'main_single_species.py'
# data and competitor
import sys
import numpy as np
import time

sys.path.append("src")
from differential_equations import lag_min, delta_max
from analytical_calculations import compute_a_and_b, compute_ap_and_bp
from simulation_functions import run_competition_in_parallel

np.random.seed(18)
constant_index = {'T0':0, 'Tab':1}
folder = 'test' #'competition_two_species'


###########################
## Simulation parameters ##
###########################
bac_res = 20                          # resolution in bacterial parameters
ab_res  = 20                          # resolution in antibiotic parameters
t_res   = 10                           # resolution in time array
tot_cycles  = 1_000



###########################
## Antibiotic parameters ##
###########################
# for chosing which time parameter to keep constant.
ic = constant_index['T0']                       # 'T0' or 'Tab'
T_const = 5                                     # value of the constant parameter

# defining parameter arrays
T_max = [12, 24]                                # upper bounds on meningful values for T0 and Tab
p_arr = np.linspace(0, 1, ab_res)               # probability of antibiotics
T_arr = np.linspace(0, T_max[1-ic], ab_res)     # time array

T_labels = ['T0', 'Tab']
T_values = [T_const, T_arr]
T0  = T_values[ic]
Tab = T_values[1-ic]



##########################
## Bacterial parameters ##
##########################
# empty arrays parameters
lag   =  np.ones([2, bac_res, bac_res]) * lag_min
delta = np.zeros([2, bac_res, bac_res])

# competitor parameters
I = np.ones(bac_res)
lag[1]   += np.outer(I, np.linspace(0, max(T0 + Tab), bac_res))
delta[1] += np.outer(np.linspace(0, delta_max, bac_res), I)

# transforming parameters to a-b space
a,  b  = compute_a_and_b(lag, delta)
ap, bp = compute_ap_and_bp(lag, delta)


########################
## Running Simulation ##
########################
tic = time.time()
bac_args = [lag, delta, a, b, ap, bp]
ab_args = [p_arr, T0, Tab]
sim_args = [ab_res, bac_res, t_res, tot_cycles]

if __name__ == '__main__':
    S_frac, lag_opt, del_opt = run_competition_in_parallel(bac_args, ab_args, sim_args)

# saving
np.savetxt(f"data/{folder}/optimal_Sfrac-{T_labels[ic]}{T_const}.txt", S_frac)
np.savetxt(f"data/{folder}/optimal_lag-{T_labels[ic]}{T_const}.txt", lag_opt)
np.savetxt(f"data/{folder}/optimal_delta-{T_labels[ic]}{T_const}.txt", del_opt)

toc = time.time()
print(f'\nTime: {(toc - tic) / 3600}h\n')
