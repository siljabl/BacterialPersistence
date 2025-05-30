# Functions for computing optimal parameters from functions in model_equations.py
import numpy as np
from differential_equations import f, S0
from analytical_calculations import solve_constants


####################
## Single species ##
####################
# consumption time of single species, as defined in thesis
# def analytical_fitness(bac_args, ab_args):
#     a = bac_args['a']
#     b = bac_args['b']
#     ap = bac_args['ap']
#     bp = bac_args['bp']

#     p   = ab_args['p']
#     T0  = ab_args['T0']
#     Tab = ab_args['Tab']

#     # consumption time without antibiotics
#     Ts = (1 / b) * np.log((a + b) / (a * b * f))    

#     prefactor = a * b * f / ((ap - bp) * (a + b) ** 2)
#     a_bp = a - bp
#     a_ap = a - ap
#     b_ap = b + ap
#     b_bp = b + bp

#     exp_aT0 = np.exp(-a * T0)
#     exp_bT0 = np.exp(b * T0)
#     exp_apT = np.exp(-ap * Tab)
#     exp_bpT = np.exp(-bp * Tab)

#     gT0 = a * b * f * (exp_bT0 - exp_aT0) / (a + b)
#     gT = prefactor * (a + b) * (a_bp * exp_bT0 + b_bp * exp_aT0) * exp_bpT - \
#          prefactor * (a + b) * (a_ap * exp_bT0 + b_ap * exp_aT0) * exp_apT

#     b_term = b_ap * (a_bp * exp_bT0 + b_bp * exp_aT0) * exp_bpT - \
#              b_bp * (a_ap * exp_bT0 + b_ap * exp_aT0) * exp_apT

#     # consumption time with antibiotics
#     Ts_ab = (1 / b) * np.log((1 + gT - gT0) / (prefactor * b_term))
    
#     # weighted average
#     Ts_avrg = (1 - p) * Ts + p * (Tab + T0 + Ts_ab)

#     return 1 / Ts_avrg



# consumption time of single species, as defined in thesis
def analytical_fitness(bac_args, ab_args):
    delta = bac_args['delta']
    a = bac_args['a']
    b = bac_args['b']

    p   = ab_args['p']
    T0  = ab_args['T0']
    Tab = ab_args['Tab']
    T   = ab_args['T']

    exp_aT0 = np.exp(-a * T0)
    exp_bT0 = np.exp(b * T0)
    exp_aT  = np.exp(-a * (T))
    exp_bT  = np.exp(b * (T))

    C1, C2 = solve_constants(bac_args, ab_args, stage="pre")
    E1, E2 = solve_constants(bac_args, ab_args, stage="post")

    int_gT0 = (C1/b)*(exp_bT0 - 1) - (C2/a)*(exp_aT0 - 1)

    # consumption time without antibiotics
    Ts    = (1 / b) * np.log((b / C1)*(S0/(1-delta) - C2/a) + 1)

    # consumption time with antibiotics
    Ts_ab = (1 / b) * np.log((b / E1)*(S0/(1-delta) - (E2/a) - int_gT0) + 1) + T
    
    # weighted average
    Ts_avrg = (1 - p) * Ts + p * Ts_ab

    return 1 / Ts_avrg




# computing optimal parameters for all combinations of antibiotic parameters
def compute_optimal_parameters(bac_args, ab_args):
    lag   = bac_args['lag']
    delta = bac_args['delta']

    p   = ab_args['p']
    T0  = ab_args['T0']
    Tab = ab_args['Tab']
    ab_res = len(p)

    # generalizing time arrays
    T0  = T0  * np.ones(ab_res)
    Tab = Tab * np.ones(ab_res)

    # output arrays
    lag_opt = np.zeros([ab_res, ab_res])
    del_opt = np.zeros([ab_res, ab_res])
    F = np.zeros([ab_res, ab_res])

    for ip in range(ab_res):  						        # probability loop
        for it in range(ab_res):  					        # duration / application time loop
            ab_arg = {'p':p[ip], 'T0':T0[it], 'Tab':Tab[it], 'T':T0[it]+Tab[it]}				# updating antibiotic parameters

            F_matrix = analytical_fitness(bac_args, ab_arg)				# computing fitness matrix corresponding to all bacterial parameters
            F_max = (F_matrix == F_matrix.max())				        # masking shortest time
            
            F[ip, it] = F_matrix[F_max]					                # saving max fitness
            lag_opt[ip, it] = lag[0, (F_max.sum(0)).astype(bool)]		# identifying corresponding lag time
            del_opt[ip, it] = delta[(F_max.sum(1)).astype(bool), 0]  	# identifying corresponding persistence

        print(100 * np.round(ip / ab_res, 4), "%")			# printing progress

    return lag_opt, del_opt, F



