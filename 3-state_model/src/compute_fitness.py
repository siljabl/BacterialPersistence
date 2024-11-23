# Functions for computing optimal parameters from functions in model_equations.py
import numpy as np

from differential_equations import S0
from estimate_consumption_time import estimate_Ts



# consumption time of single species, as defined in thesis
def analytical_fitness(eq_params, ab_params):
    '''
    Computing fitness, weighted by probability of antiobiotic cycle.
    '''
    #b  = eq_params['b']
    p  = ab_params['p']
    T  = ab_params['T']
    #T0 = ab_params['T0']
    
    #B,_,_ = solve_constants(eq_params, ab_params, 'post')

    Ts_min = 0 # np.log(1/f)
    Ts_max = 100 #np.max(T0 * 40)
    #Ts_arr = np.linspace(Ts_min, Ts_max, 100)

    Ts, Ts_p = estimate_Ts(eq_params, ab_params, iter=3)

    # weighted average
    Ts_avrg = (1 - p) * Ts + p * (T + Ts_p)

    return 1 / Ts_avrg



def transform_fitness_to_bac_parameters(fitness, bac_params):
    '''
    This function identifies the optimal parameters from the maximal fitness
    '''
    λd = bac_params['λd']
    λr = bac_params['λr']
    δ  = bac_params['δ']

    F_max = np.max(fitness)
    idx = (fitness == F_max)

    return F_max, λd[idx], λr[idx], δ[idx]
