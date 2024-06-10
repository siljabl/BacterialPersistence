import sys
import numpy as np
import matplotlib.pyplot as plt
import time

sys.path.append('src')
from initialization import initialize_bacterial_parameter_arrays, initialise_system
from evolution import evolve_system
from test_plots import *
#from model_equations import initialise_system, evolve_system


###########################
## Simulation parameters ##
###########################
# np.random.seed(2)
sim_params = {'mutation':   False, 
              'extinct':    False,
              'mut_seed':   'min', 
              'tot_cycles': 10_000}


###########################
## Antibiotic parameters ##
###########################
T0  = 0.0
Tab = 12.0
T   = T0 + Tab
T0_max = 12

p_arr = np.linspace(0.1, 0.9, 5)
p = 0.1

ab_params = {'p':p, 'T0':T0, 'Tab':Tab, 'T':T0+Tab, 'T0_max':T0_max}

dλ = 1                  	# wake up rate from dormancy
dδ = 0.001               	# rate of spontaneous persistence
Ɛ = 10**(-3)

################
## Simulation ##
################
r_arr = np.random.rand(sim_params['tot_cycles'])
sim_params['r_arr'] = r_arr

bac_params = initialize_bacterial_parameter_arrays(ab_params, dλ, dδ)
bac_params['Ɛ'] = Ɛ
plot_bacterial_parameters(bac_params)

populations = initialise_system(bac_params, sim_params)
λd, λr, δ, p_dominant, n_extinct, p_dists, cycles = evolve_system(populations, bac_params, ab_params, sim_params)
plot_cycles(cycles, bac_params, sim_params)

np.savetxt(f"data/competition_average/average_λd-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt", λd)
np.savetxt(f"data/competition_average/average_λr-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt", λr)
np.savetxt(f"data/competition_average/average_δ-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt",  δ)
np.savetxt(f"data/competition_average/dominant_species-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt", p_dominant)
np.savetxt(f"data/competition_average/number_of_extinctions-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt", n_extinct)
np.savetxt(f"data/competition_average/population_distributions-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt", p_dists) 
