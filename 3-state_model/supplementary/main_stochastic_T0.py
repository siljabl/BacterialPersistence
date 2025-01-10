import sys
import argparse
import numpy as np

sys.path.append("../src")
from differential_equations  import λ_min, δ_max
from analytical_calculations import compute_a_and_b, compute_ap_and_bp
from compute_fitness         import analytical_fitness, transform_fitness_to_bac_parameters


parser = argparse.ArgumentParser(description='Run and plot fitness only for one set of antibiotic parameters.')
parser.add_argument('p',   type=float, help='probability of antibiotics')
parser.add_argument('T0',  type=float, help='application time of antibiotics')
parser.add_argument('Tab', type=float, help='duration of antibiotics')
parser.add_argument('-bac_res', type=int, help='resolution on bacterial parameters',  nargs='?', default=100)
args = parser.parse_args()

p_spon = np.linspace(0, 1, 11)


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
ab_params_trig  = {'p':args.p, 'T0':0, 'Tab':args.Tab, 'T':args.Tab}
ab_params_spon  = {'p':args.p, 'T0':args.T0, 'Tab':args.Tab, 'T':args.T0+args.Tab}

###########################################
## Looping through antibiotic parameters ##
###########################################
F  = np.zeros(len(p_spon))
δ  = np.zeros(len(p_spon))
λd = np.zeros(len(p_spon))
λr = np.zeros(len(p_spon))


fitness_trig = analytical_fitness(eq_params, ab_params_trig)
fitness_spon = analytical_fitness(eq_params, ab_params_spon)

for i in range(len(p_spon)):
    fitness = (1-p_spon[i])*fitness_trig + p_spon[i]*fitness_spon
    optimal_params = transform_fitness_to_bac_parameters(fitness, bac_params)

    F[i]  = optimal_params[0]
    λd[i] = optimal_params[1][0]
    λr[i] = optimal_params[2][0]
    δ[i]  = optimal_params[3][0]

    if optimal_params[3][0] == 0:
        λr[i] = 0

    print(i)


np.savetxt(f"../data/supplementary/F_stochasticT0_p{args.p}_T0max{args.T0}_Tab{args.Tab}.txt", F)
np.savetxt(f"../data/supplementary/δ_stochasticT0_p{args.p}_T0max{args.T0}_Tab{args.Tab}.txt",  δ)
np.savetxt(f"../data/supplementary/λd_stochasticT0_p{args.p}_T0max{args.T0}_Tab{args.Tab}.txt", λd)
np.savetxt(f"../data/supplementary/λr_stochasticT0_p{args.p}_T0max{args.T0}_Tab{args.Tab}.txt", λr)
