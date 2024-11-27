# computing optimal parameter of single species with restricted nutrients
import sys
import argparse
import numpy as np
from datetime import datetime

sys.path.append("src")
from config_functions        import save_config
from differential_equations  import λ_min, δ_max
from analytical_calculations import compute_a_and_b, compute_ap_and_bp
from compute_fitness         import analytical_fitness, transform_fitness_to_bac_parameters


#####################
## Input arguments ##
#####################
parser = argparse.ArgumentParser(description='NULL')
parser.add_argument('folder',          type=str, help="Folder for saving data.")
parser.add_argument('T0',              type=int, help='Time at which antibiotics are applied.')
parser.add_argument('Tab_max',         type=int, help='Upper limit on antibiotic duration.')
parser.add_argument('-Tab_min',        type=int, help='Lower limit on antibiotic duration.', nargs='?', default=0)
parser.add_argument('-antibiotic_res', type=int, help='resolution on antibiotic parameters', nargs='?', default=41)
parser.add_argument('-bacterial_res',  type=int, help='resolution on bacterial parameters',  nargs='?', default=100)
args = parser.parse_args()

folder  = args.folder
T0      = args.T0
Tab_max = args.Tab_max
Tab_min = args.Tab_min
ab_res  = args.antibiotic_res
bac_res = args.bacterial_res
stoch_param = 'none'
Tab_res = 11


###########################
## Antibiotic parameters ##
###########################
p_arr   = np.linspace(0, 1, ab_res)
Tab_arr = np.linspace(0, Tab_max, Tab_res)


##########################
## Bacterial parameters ##
##########################
λ_arr = np.linspace(0, T0 + Tab_max + 1, bac_res) + λ_min
δ_arr = np.linspace(0, δ_max, bac_res)

λd = np.outer(np.ones(bac_res), np.outer(np.ones(bac_res), λ_arr)).reshape(bac_res, bac_res, bac_res) - λ_min / 1000    # avoid overflow by distinguishing λd and λr
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
λd_opt = np.zeros([ab_res, Tab_res])
λr_opt = np.zeros([ab_res, Tab_res])
δ_opt  = np.zeros([ab_res, Tab_res])
F_max  = np.zeros([ab_res, Tab_res])


i = 0
for p in p_arr:
    j = 0
    for Tab in Tab_arr:
        ab_params = {'p':p, 'T0':T0, 'Tab':Tab, 'T':T0+Tab}
        fitness = analytical_fitness(eq_params, ab_params)
        optimal_params = transform_fitness_to_bac_parameters(fitness, bac_params)
        
        F_max[i,j]  = optimal_params[0]
        λd_opt[i,j] = optimal_params[1][0]
        λr_opt[i,j] = optimal_params[2][0]
        δ_opt[i,j]  = optimal_params[3][0]

        if optimal_params[3][0] == 0:
            λr_opt[i,j] = 0

        j = j + 1
    
    print(p)
    i = i + 1


np.savetxt(f"{folder}/single_optimal_F_max.txt", F_max)
np.savetxt(f"{folder}/single_optimal_λd-T0_{T0}.txt", λd_opt)
np.savetxt(f"{folder}/single_optimal_λr-T0_{T0}.txt", λr_opt)
np.savetxt(f"{folder}/single_optimal_δ-T0_{T0}.txt", δ_opt)


config = {"date"    : datetime.today().strftime('%Y-%m-%d'),
          "ab_res"  : ab_res,
          "bac_res" : bac_res,
          "T0"      : T0,
          "Tab_min" : Tab_min,
          "Tab_max" : Tab_max,
          "λd_min"  : np.min(λd),
          "λd_max"  : np.max(λd),
          "λr_min"  : np.min(λr),
          "λr_max"  : np.max(λr),
          "δ_min"   : np.min(δ),
          "δ_max"   : np.max(δ)}

save_config(config, folder)
