import sys
import argparse
import numpy as np

sys.path.append("../src")
from differential_equations  import λ_min, δ_max
from analytical_calculations import compute_a_and_b, compute_ap_and_bp
from compute_fitness         import analytical_fitness


parser = argparse.ArgumentParser(description='Run and plot fitness only for one set of antibiotic parameters.')
parser.add_argument('T0',  type=float, help='application time of antibiotics')
parser.add_argument('Tab', type=float, help='duration of antibiotics')
parser.add_argument('-bac_res', type=int, help='resolution on bacterial parameters',  nargs='?', default=100)
args = parser.parse_args()


p_arr = np.linspace(0.1, 0.9, 9)

##########################
## Bacterial parameters ##
##########################
λ_arr = np.linspace(0, args.T0 + args.Tab + 1, args.bac_res) + λ_min
δ_arr = np.linspace(0, δ_max, args.bac_res)

λd = np.outer(np.ones(args.bac_res), np.outer(np.ones(args.bac_res), λ_arr)).reshape(args.bac_res, args.bac_res, args.bac_res) - λ_min / 1000    # avoid overflow by distinguishing λd and λr
λr = np.outer(np.outer(np.ones(args.bac_res), λ_arr), np.ones(args.bac_res)).reshape(args.bac_res, args.bac_res, args.bac_res)
δ  = np.outer(np.outer(δ_arr, np.ones(args.bac_res)), np.ones(args.bac_res)).reshape(args.bac_res, args.bac_res, args.bac_res)


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
F_trig = np.zeros(len(p_arr))
F_spon = np.zeros(len(p_arr))


for i in range(len(p_arr)):
    ab_params  = {'p':p_arr[i], 'T0':args.T0, 'Tab':args.Tab, 'T':args.T0+args.Tab}
    fitness = analytical_fitness(eq_params, ab_params)

    F_trig[i] = np.max(fitness[0,0,:])
    F_spon[i] = np.max(fitness[:,:,0])

    print(i / (len(p_arr)-1))

    
F = np.array([F_trig, F_spon]).T
np.savetxt(f"../data/supplementary/strategy_fitness_T0{args.T0}_Tab{args.Tab}.txt", F, header="triggered spontaneous")
