import sys
import pickle
import argparse
import numpy as np

sys.path.append('src')
from initialization import initialize_bacterial_parameter_arrays, initialise_system
from test_plots import *

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


##########################
## Bacterial parameters ##
##########################
dλ = 1                  	# wake up rate from dormancy
dδ = 0.001               	# rate of spontaneous persistence
Ɛ = 10**(-3)

bac_params = initialize_bacterial_parameter_arrays(ab_params, dλ, dδ)
bac_params['Ɛ'] = Ɛ


#####################
## Plot parameters ##
#####################
plot_bacterial_parameters(bac_params, folder)

#################
## Plot cycles ##
#################
with open(f"{folder}/solve_cycles.pkl", 'rb') as file:
    cycles = pickle.load(file)

plot_cycles(cycles, bac_params, ab_params, sim_params, folder)