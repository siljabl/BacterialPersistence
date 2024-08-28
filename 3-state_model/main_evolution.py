import sys
import pickle
import argparse
import numpy as np
from datetime import datetime

sys.path.append('src')
from initialization import initialize_bacterial_parameter_arrays, initialise_system
from evolution import evolve_system
# from test_plots import *

parser = argparse.ArgumentParser(description='Competition between N species for tot_cycles cycles.')
parser.add_argument('folder',          type=str, help="Folder for saving data.")
parser.add_argument('p',   type=float, help='frequency of antibiotics.')
parser.add_argument('T0',  type=float, help='application time of antibiotics')
parser.add_argument('Tab', type=float, help='duration of antibiotics')
parser.add_argument('--tot_cycles', type=int,  nargs='?', \
                    help='number of feast-famine cycles that are performed', default=10_000)
args = parser.parse_args()


###########################
## Simulation parameters ##
###########################
# np.random.seed(2)
tot_cycles = args.tot_cycles
sim_params = {'mutation':   False, 
              'extinct':    False,
              'mut_seed':   'min', 
              'tot_cycles': tot_cycles}


###########################
## Antibiotic parameters ##
###########################
folder = args.folder
p      = args.p
T0     = args.T0
Tab    = args.Tab
T      = T0 + Tab
T0_max = 12


ab_params = {'p':p, 'T0':T0, 'Tab':Tab, 'T':T0+Tab, 'T0_max':T0_max}

dλ = 1                  	# wake up rate from dormancy
dδ = 0.001               	# rate of spontaneous persistence
Ɛ = 10**(-3)

################
## Simulation ##
################
r_arr = np.random.rand(sim_params['tot_cycles'])
#r_arr[1000:] = 0.5
sim_params['r_arr'] = r_arr

bac_params = initialize_bacterial_parameter_arrays(ab_params, dλ, dδ)
bac_params['Ɛ'] = Ɛ


populations = initialise_system(bac_params, sim_params)
λd, λr, δ, p_dominant, n_extinct, p_dists, cycles = evolve_system(populations, bac_params, ab_params, sim_params)


np.savetxt(f"{folder}/competition_average_λd-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt", λd)
np.savetxt(f"{folder}/competition_average_λr-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt", λr)
np.savetxt(f"{folder}/competition_average_δ-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt",  δ)
np.savetxt(f"{folder}/r_arr_specified_random_array-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt",           r_arr)

with open(f"{folder}/solve_cycles_p_{p:0.1f}.pkl", "wb") as file: 
    pickle.dump(cycles, file) 

# np.savetxt(f"data/competition_average/dominant_species-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt", p_dominant)
# np.savetxt(f"data/competition_average/number_of_extinctions-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt", n_extinct)
# np.savetxt(f"data/competition_average/population_distributions-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt", p_dists) 
