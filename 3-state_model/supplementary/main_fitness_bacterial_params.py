import sys
import argparse
import numpy as np

sys.path.append("../src")
from differential_equations  import λ_min, δ_max
from analytical_calculations import compute_a_and_b, compute_ap_and_bp
from compute_fitness         import analytical_fitness


parser = argparse.ArgumentParser(description='Run and plot fitness only for one set of antibiotic parameters.')
parser.add_argument('p',   type=float, help='probability of antibiotics')
parser.add_argument('T0',  type=float, help='application time of antibiotics')
parser.add_argument('Tab', type=float, help='duration of antibiotics')
parser.add_argument('-bac_res', type=int, help='resolution on bacterial parameters',  nargs='?', default=100)
args = parser.parse_args()



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
ab_params  = {'p':args.p, 'T0':args.T0, 'Tab':args.Tab, 'T':args.T0+args.Tab}

###########################################
## Looping through antibiotic parameters ##
###########################################
F  = np.zeros([args.bac_res, args.bac_res])
δ  = np.zeros([args.bac_res, args.bac_res])
λd = np.zeros([args.bac_res, args.bac_res])
λr = np.zeros([args.bac_res, args.bac_res])

fitness = analytical_fitness(eq_params, ab_params)
#fitness_projection = np.max(fitness, axis=0)

for i in range(args.bac_res):
    for j in range(args.bac_res):
        idx = np.argmax(fitness[:,i,j])

        F[i,j]  = fitness[idx,i,j]
        δ[i,j]  = bac_params['δ'][idx,i,j]
        λd[i,j] = bac_params['λd'][idx,i,j]
        λr[i,j] = bac_params['λr'][idx,i,j]



#idx = (fitness == fitness_projection)


np.savetxt(f"../data/supplementary/fitness_δprojection_p{args.p}_T0{args.T0}_Tab{args.Tab}.txt", F)
np.savetxt(f"../data/supplementary/δ_δprojection_p{args.p}_T0{args.T0}_Tab{args.Tab}.txt",  δ)
np.savetxt(f"../data/supplementary/λd_δprojection_p{args.p}_T0{args.T0}_Tab{args.Tab}.txt", λd)
np.savetxt(f"../data/supplementary/λr_δprojection_p{args.p}_T0{args.T0}_Tab{args.Tab}.txt", λr)
