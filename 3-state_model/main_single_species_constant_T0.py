# computing optimal parameter of single species with restricted nutrients
import sys
import time
import argparse
import numpy as np

sys.path.append("src")
from differential_equations  import p_max, λ_min, δ_max, S0
from analytical_calculations import compute_a_and_b, compute_ap_and_bp
from compute_fitness import analytical_fitness, transform_fitness_to_bac_parameters

# Plotting configurations
# font = {'family': 'Times New Roman',
#         'weight': 'normal',
#         'size': 20}
# mpl.rc('font', **font)

#λ_cmap = mpl.cm.get_cmap('viridis')
#δ_cmap = mpl.cm.get_cmap('plasma')

#####################
## Input arguments ##
#####################
parser = argparse.ArgumentParser(description='NULL')
parser.add_argument('T0',              type=int, help='Time at which antibiotics are applied.')
parser.add_argument('Tab_max',         type=int, help='Upper limit on antibiotic duration.')
parser.add_argument('-antibiotic_res', type=int, help='resolution on antibiotic parameters', nargs='?', default=11)
parser.add_argument('-bacterial_res',  type=int, help='resolution on bacterial parameters',  nargs='?', default=10)
args = parser.parse_args()

T0      = args.T0
Tab_max = args.Tab_max
ab_res  = args.antibiotic_res
bac_res = args.bacterial_res
stoch_param = 'none'


###########################
## Antibiotic parameters ##
###########################
p_arr   = np.linspace(0, 1, ab_res)       # probability of antibiotics
Tab_arr = np.linspace(0, Tab_max, ab_res)     # time array


##########################
## Bacterial parameters ##
##########################
λ_arr = np.linspace(0, T0 + Tab_max, bac_res) + λ_min
δ_arr = np.linspace(0, δ_max, bac_res)

λd = np.outer(np.ones(bac_res), np.outer(np.ones(bac_res), λ_arr)).reshape(bac_res, bac_res, bac_res) - λ_min / 1000
λr = np.outer(np.outer(np.ones(bac_res), λ_arr), np.ones(bac_res)).reshape(bac_res, bac_res, bac_res)
δ  = np.outer(np.outer(δ_arr, np.ones(bac_res)), np.ones(bac_res)).reshape(bac_res, bac_res, bac_res)


###########################
## Simulation Parameters ##
###########################
a, b   = compute_a_and_b(λr, δ)
ap, bp = compute_ap_and_bp(λr, δ)
c = 1/λd

bac_params = {'λd':λd, 'λr':λr, 'δ':δ}
eq_params  = {'a':a, 'b':b, 'c':c, 'ap':ap, 'bp':bp}


###########################################
## Looping through antibiotic parameters ##
###########################################
λd_opt = np.zeros([ab_res, ab_res])
λr_opt = np.zeros([ab_res, ab_res])
δ_opt  = np.zeros([ab_res, ab_res])
F_max  = np.zeros([ab_res, ab_res])


i = 0
for p in p_arr:
    j = 0
    for Tab in Tab_arr:
        ab_params = {'p':p, 'T0':T0, 'Tab':Tab, 'T':T0+Tab}
        fitness = analytical_fitness(eq_params, ab_params)
        optimal_params = transform_fitness_to_bac_parameters(fitness, bac_params)
        
        F_max[i,j] = optimal_params[0]
        λd_opt[i,j] = optimal_params[1][0]
        λr_opt[i,j] = optimal_params[2][0]
        δ_opt[i,j]  = optimal_params[3][0]

        j = j + 1
    i = i + 1



np.savetxt(f"data/optimal_λd-T0_{T0}.txt", λd_opt)
np.savetxt(f"data/optimal_λr-T0_{T0}.txt", λr_opt)
np.savetxt(f"data/optimal_δ-T0_{T0}.txt", δ_opt)
