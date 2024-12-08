# Functions for computing optimal parameters from functions in model_equations.py
import numpy as np
import psutil
from multiprocessing.pool import Pool

from optimal_from_file   import optimal_parameters_from_data
from competition_fitness import run_competition



# Looping through antibiotic parameters
def looping_through_antibiotic_parameters(bac_args, ab_args, sim_args):
    p_arr, T0, Tab = ab_args
    ab_res, save_data = sim_args[0], sim_args[5]
    #np.random.seed(int(min(p_arr) * ab_res))            # setting different random seeds for each subdomain (parallel simulation)
    
    T0 = T0 * np.ones(ab_res)                           # preparing time array
    Tab = Tab * np.ones(ab_res)                         # preparing time array
    opt_params = np.zeros([len(p_arr), ab_res, 3])  	# output array

    ip = 0
    # Probability loop
    for p in p_arr:
        # Time loop
        for it in range(ab_res):
            ab_args = [p, T0[it], Tab[it]]        # subset of antibiotic parameters
            bac_args = optimal_parameters_from_data(bac_args, ab_args)     # compute winner parameters

            # Without extinction
            opt_params[ip, it], prob_ext = run_competition(bac_args, ab_args, sim_args)[0:2]
            if save_data:
                config = str(int(10*p)*10) + "-T0" + str(int(T0[it])) + "-Tab" + str(int(Tab[it]))
                np.savetxt("data/extinction_frequency/optimal_extinction_prob-p" + config, prob_ext[0])
                np.savetxt("data/extinction_frequency/competitor_extinction_prob-p" + config, prob_ext[1])

            if it % 10 == 0 and 100 * p % 10 == 0:
                print(100 * np.round((it + 10) / ab_res, 2), "% of p = " + str(np.round(p, 2)))  # print progression
        ip += 1

    return opt_params, p_arr.min()



# run simulation
def run_competition_in_parallel(bac_args, ab_args, sim_args):
    p_arr, T0, Tab = ab_args
    ab_res, bac_res = sim_args[0:2]

    cores = psutil.cpu_count(logical=False)             # number of available cores
    print(f'Available cores: {cores}')
    # Domain decomposition
    width = int(ab_res / cores)                         # width of subdomain
    domain_order = np.zeros(cores)                      # array for sorting domains
    results = np.zeros([width, ab_res, 3, cores])       # output array

    ab_subdomains = []                                  # array of subdomains
    for i in range(cores):
        ab_subdomains.append([p_arr[i * width: (i + 1) * width], T0, Tab])

    # Running jobs in parallel
    with Pool(cores) as pool:
        jobs = [pool.apply_async(looping_through_antibiotic_parameters, (bac_args, ab_subdomain, sim_args)) for ab_subdomain in ab_subdomains]
        ir = 0
        # Collecting results
        for result in [job.get() for job in jobs]:
            results[:, :, :, ir], domain_order[ir] = result
            ir += 1

    # Assemble domain
    idx_sorted = np.argsort(domain_order)
    S_opt   = np.concatenate([results[:, :, 0, i] for i in idx_sorted])
    lag_opt = np.concatenate([results[:, :, 1, i] for i in idx_sorted])
    del_opt = np.concatenate([results[:, :, 2, i] for i in idx_sorted])

    return S_opt, lag_opt, del_opt
