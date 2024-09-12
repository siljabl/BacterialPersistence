# Functions for computing optimal parameters from functions in model_equations.py
import numpy as np
import psutil
from multiprocessing.pool import Pool

from model_equations import a_b, ap_bp
from model_equations import analytical_population, analytical_growth, analytical_decay
from model_equations import n0, S0, f, n_min

import matplotlib.pyplot as plt


####################
## Single species ##
####################
# consumption time of single species, as defined in thesis
def analytical_fitness(bac_args, ab_args):
    _, _, a, b, ap, bp = bac_args
    p, T0, Tab = ab_args

    # consumption time without antibiotics
    Ts = (1 / b) * np.log((a + b) / (a * b * f))    

    prefactor = a * b * f / ((ap - bp) * (a + b) ** 2)
    a_bp = a - bp
    a_ap = a - ap
    b_ap = b + ap
    b_bp = b + bp

    exp_aT0 = np.exp(-a * T0)
    exp_bT0 = np.exp(b * T0)
    exp_apT = np.exp(-ap * Tab)
    exp_bpT = np.exp(-bp * Tab)

    gT0 = a * b * f * (exp_bT0 - exp_aT0) / (a + b)
    gT = prefactor * (a + b) * (a_bp * exp_bT0 + b_bp * exp_aT0) * exp_bpT - \
         prefactor * (a + b) * (a_ap * exp_bT0 + b_ap * exp_aT0) * exp_apT

    b_term = b_ap * (a_bp * exp_bT0 + b_bp * exp_aT0) * exp_bpT - \
             b_bp * (a_ap * exp_bT0 + b_ap * exp_aT0) * exp_apT

    # consumption time with antibiotics
    Ts_ab = (1 / b) * np.log((1 + gT - gT0) / (prefactor * b_term))
    
    # weighted average
    Ts_avrg = (1 - p) * Ts + p * (Tab + T0 + Ts_ab)

    return 1 / Ts_avrg


# computing optimal parameters for all combinations of antibiotic parameters
def compute_optimal_parameters(bac_args, ab_args):
    lag, delta = bac_args[0:2]
    p, T0, Tab = ab_args
    ab_res, bac_res = len(p), len(lag)

    # generalizing time arrays
    T0 = T0 * np.ones(ab_res)
    Tab = Tab * np.ones(ab_res)

    # output arrays
    lag_opt = np.zeros([ab_res, ab_res])
    del_opt = np.zeros([ab_res, ab_res])
    F = np.zeros([ab_res, ab_res])

    for ip in range(ab_res):  						        # probability loop
        for it in range(ab_res):  					        # duration / application time loop
            ab_arg = [p[ip], T0[it], Tab[it]]				# updating antibiotic parameters

            F_matrix = analytical_fitness(bac_args, ab_arg)				# computing fitness matrix corresponding to all bacterial parameters
            F_max = (F_matrix == F_matrix.max())				        # masking shortest time
            #np.savetxt("../../../../data2", F_matrix)

            F[ip, it] = F_matrix[F_max]					                # saving max fitness
            lag_opt[ip, it] = lag[0, (F_max.sum(0)).astype(bool)]		# identifying corresponding lag time
            del_opt[ip, it] = delta[(F_max.sum(1)).astype(bool), 0]  	# identifying corresponding persistence

        print(100 * np.round(ip / ab_res, 4), "%")			# printing progress

    return lag_opt, del_opt, F_matrix



#################
## Competition ##
#################
# Function for finding ts = t: S(t) = 0
def Ts_approximation(t, n_t0, a, b):
    d_t0, g_t0, S_0 = n_t0                                          # initial values
    form = np.shape(S_0)                                            # shape of output
    NA = np.newaxis                                                 # extra axis for array multiplication

    t_steps = len(t)                                                # length of time array
    dt_new = t[1] - t[0]                                            # time step
    t_rand = (np.random.rand(form[0], form[1]) - 0.5) * dt_new  	# random offset

    # increasing precision on Ts
    for i in range(3):
        t_temp = t[:, NA, NA, NA] + t_rand                          # temporary time array
        exp_bt = np.exp(b * t_temp)                                 # time dependent terms
        exp_at = np.exp(-a * t_temp)                                # time dependent terms
        
        # analytical population
        g_t = d_t0 * a * b * (exp_bt - exp_at) / (a + b) + \
              g_t0 * (b * (1 + a) * exp_bt + a * (1 - b) * exp_at) / (a + b)  	# size of growing population
        S = S_0 - (g_t - g_t0).sum(axis=len(np.shape(g_t)) - 2 - 1)             # amount of substrate
        #print("S:", np.shape(S))

        # Finding upper and lower limit on ts
        idx = np.where(np.diff(np.sign(S), axis=0) < 0)             # identify indices where S changes sign
        #print("diff:", np.diff(np.sign(S), axis=0).shape)
        #print("idx:", np.shape(idx))
        idx_sorted = np.lexsort((idx[2], idx[1]))                   # sorting indices
        #print("idx_sorted:", np.shape(idx_sorted))

        Ts_lowlim = t[idx[0]]                                        # indentifying corresponding time
        Ts_lowlim = Ts_lowlim[idx_sorted].reshape(form) + t_rand     # adding offset
        t_rand = Ts_lowlim                                           # updating offset

        # Preparing for next round
        t = np.linspace(0, dt_new, t_steps)                         # new time array
        dt_old = dt_new                                             # saving time step
        dt_new = t[1] - t[0]                                        # updating time step

    Ts_above = Ts_lowlim                                            # corresponding time
    S_above = S[idx]                                                # corresponding substrate
    S_below = S[idx[0] + 1, idx[1], idx[2]]                         # value of substrate after sign change
    S_above = S_above[idx_sorted].reshape(form)                     # reshaping substrate to form of output
    S_below = S_below[idx_sorted].reshape(form)                     # reshaping substrate to form of output

    # Linear approximation around S = 0
    dS = S_below - S_above                                          # change in substrate
    slope = dS / dt_old                                             # computing slope
    const = S_above - slope * Ts_above                              # computing offset
    Ts = - const / slope                                            # ts: S(ts) = 0

    return Ts



# Function that yields optimal parameters of single species winner
def optimal_parameters_from_data(bac_args, ab_args):

    lag, delta, a, b, ap, bp = bac_args
    p, T0, Tab = ab_args
    ab_res = len(np.loadtxt("data/low_resolution/optimal_lag-T0" + str(int(T0))))
    ip = int(p * (ab_res-1))
    it = int(Tab * (ab_res-1) / 24)

    lag[0] = np.loadtxt("data/low_resolution/optimal_lag-T0" + str(int(T0)))[ip, it]
    delta[0] = np.loadtxt("data/low_resolution/optimal_delta-T0" + str(int(T0)))[ip, it]

    # Transforming to a-b scheme
    a[0], b[0] = a_b(lag[0], delta[0])
    ap[0], bp[0] = ap_bp(lag[0], delta[0])

    return [lag, delta, a, b, ap, bp]


# Finding optimal competitor parameters
def run_competition(bac_args, ab_args, sim_args):
    # Sorting inputs
    lag, delta, a, b, ap, bp = bac_args
    p, T0, Tab = ab_args[0:3]
    _, bac_res, t_res, tot_cycles, reps, _ = sim_args

    t_min, t_max = 12 - T0, 15 + (T0 + Tab)                     # time limits
    t_arr = np.linspace(t_min, t_max, t_res)                    # time array

    S_frac_cycle = np.zeros([tot_cycles, 2, bac_res, bac_res])  # array for output
    p_extinct = np.zeros_like(lag)  # array for output

    # Compute ensemble average
    for rep in range(reps):
        # Initial populations, as [d(0), g(0), S(0)]
        n_0 = [f * S0 * np.ones_like(lag), np.zeros_like(lag), S0 * np.ones_like(lag)[0]]
        extinct = np.zeros([bac_res, bac_res])                  # array for counting first extinctions
        ext = np.ones_like(lag)                                 # masking 
        r_arr = np.random.rand(tot_cycles)                      # random array for ab

        for ic in range(tot_cycles):

            # With antibiotics
            if r_arr[ic] < p:
                n_T0 = analytical_growth(T0, n_0, a, b)             # population before AB
                n_T = analytical_decay(Tab, n_T0, a, b, ap, bp)     # population after AB
                d_dead = n_T[0] < n_min                             # checking if dormant species is killed by AB
                g_dead = n_T[1] < n_min                             # checking if growing species is killed by AB
                p_extinct[d_dead * g_dead * (extinct == 0)] += 1    # counting number of species that went extinct
                ext[d_dead * g_dead * (extinct == 0)] = 0
                for i in range(2):
                    extinct[d_dead[i] * g_dead[i]] = 1              # updating counter of first extinction

                t = Ts_approximation(t_arr, n_T, a, b)              # computing ts: S(ts) = 0
                n_t = analytical_growth(t, n_T, a, b)               # population at ts
                # S_frac_cycle[ic] += (n_t[1] - n_T[1] + n_T0[1]) / ((S0 - n_t[2]) * reps)
                S_frac_cycle[ic] += (n_t[1] - n_T[1] + n_T0[1]) * ext / (S0 * reps) + (S0 - n_t[1] + n_T[1] - n_T0[1]) * ext * extinct / (S0 * reps)

            # Without antibiotics
            else:
                t = Ts_approximation(t_arr, n_0, a, b)              # computing ts: S(ts) = 0
                n_t = analytical_growth(t, n_0, a, b)               # population at ts
                # S_frac_cycle[ic] += n_t[1] / ((S0 - n_t[2]) * reps)
                S_frac_cycle[ic] += n_t[1] * ext / (S0 * reps) + (S0 - n_t[1]) * ext * extinct / (S0 * reps)

            # Preparing for next cycle
            d_0 = f * (n_t[0] + n_t[1])                             # enter dormancy and dilute
            dead = d_0 < n_min                                      # checking if species killed by dilution
            p_extinct[dead * (extinct == 0)] += 1                   # counting extinction
            ext[dead * (extinct == 0)] = 0
            for i in range(2):
                extinct[dead[i]] = 1                                # updating counter of first extinction
            n_0 = [d_0 * np.ones_like(lag), np.zeros_like(lag), S0 * np.ones_like(lag)[0]]

    # Finding optimal set of bacterial parameters
    S_frac_mean = S_frac_cycle[:, 1].mean(axis=0)                   # taking cycle average of consumption fraction
    S_max = (S_frac_mean == S_frac_mean.max())                      # finding max consumption fraction

    S_frac = S_frac_mean.max()  # saving max consumption fraction
    lag_opt = lag[1, 0][(S_max.sum(0)).astype(bool)][0]             # saving corresponding lag time
    del_opt = delta[1][(S_max.sum(1)).astype(bool), 0][0]           # saving corresponding type-II fraction

    p_extinct /= reps                                               # computing fraction of extinctions

    return np.array([S_frac, lag_opt, del_opt]), p_extinct, n_t, t, S_frac_cycle, r_arr


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

    # Verifying that domain decomposition works
    # fig, ax = plt.subplots(cores, 1, figsize=(12, 12), sharex=True)
    # for i in range(cores):
    #         ax[cores - 1 - i].imshow(results[:, :, 1, i], vmax=np.max(results[:,:,1]), origin='lower')
    # fig.show()
    # fig.savefig("figs/" + path + "/decomposition_check.png")

    return S_opt, lag_opt, del_opt
