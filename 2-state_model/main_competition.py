# Computing consumption fractions from competition between species with optimal parameters from 'main_single_species.py'
# data and competitor
import sys
import numpy as np
<<<<<<< HEAD
#import matplotlib as mpl
#import matplotlib.pyplot as plt
=======
>>>>>>> ad7101d (prepared main_competition to run on server)
import time

sys.path.append("src")
from model_equations import a_b, ap_bp, lag_min, delta_max
from simulation_functions import run_competition_in_parallel

np.random.seed(18)
<<<<<<< HEAD

#font = {'family': 'Times New Roman',
#        'weight': 'normal',
#        'size': 20}
#mpl.rc('font', **font)

#lag_cmap = mpl.cm.get_cmap('viridis')
#del_cmap = mpl.cm.get_cmap('plasma')
=======
>>>>>>> ad7101d (prepared main_competition to run on server)
constant_index = {'T0':0, 'Tab':1}



###########################
## Simulation parameters ##
###########################
bac_res = 100                          # resolution in bacterial parameters
ab_res  = 100                           # resolution in antibiotic parameters
t_res   = 10                             # resolution in time array
tot_cycles  = 10_000
repetitions = 1                        # number of repetitions for ensemble average

# data = 'new'                                # 'new' - generate and plot new data. 'old' plot old data
# save_fig = False
# save_data = False



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
lag   = np.ones([2, bac_res, bac_res]) * lag_min
delta = np.zeros([2, bac_res, bac_res])

# competitor parameters
I = np.ones(bac_res)
lag[1]   += np.outer(I, np.linspace(0, max(T0 + Tab), bac_res))
delta[1] += np.outer(np.linspace(0, delta_max, bac_res), I)

# transforming parameters to a-b space
a,  b  = a_b(lag, delta)
ap, bp = ap_bp(lag, delta)



########################
## Running Simulation ##
########################
tic = time.time()
bac_args = [lag, delta, a, b, ap, bp]
ab_args = [p_arr, T0, Tab]
sim_args = [ab_res, bac_res, t_res, tot_cycles, repetitions, False]

if __name__ == '__main__':
    S_frac, lag_opt, del_opt = run_competition_in_parallel(bac_args, ab_args, sim_args)

# saving
np.savetxt(f"data/competition_two_species/competition_Sfrac-{T_labels[ic]}{T_const}.txt", S_frac)
np.savetxt(f"data/competition_two_species/competition_lag-{T_labels[ic]}{T_const}.txt", lag_opt)
np.savetxt(f"data/competition_two_species/competition_delta-{T_labels[ic]}{T_const}.txt", del_opt)

toc = time.time()
print(f'\nTime: {(toc - tic) / 3600}h\n')
