# Functions for computing optimal parameters from functions in model_equations.py
import numpy as np

from differential_equations import S0
from estimate_consumption_time import estimate_Ts
from analytical_calculations import solve_constants


# consumption time of single species, as defined in thesis
def analytical_fitness(bac_params, eq_params, ab_params):
    '''
    Computing fitness, weighted by probability of antiobiotic cycle.
    '''
    #delta = eq_params['delta']
    a = eq_params['a']
    b = eq_params['b']
    c = eq_params['c']
    δ = bac_params['δ']

    p   = ab_params['p']
    T0  = ab_params['T0']
    Tab = ab_params['Tab']
    T   = ab_params['T']

    exp_aT0 = np.exp(-a * T0)
    exp_bT0 = np.exp(b  * T0)
    exp_cT0 = np.exp(-c * T0)
    exp_aT  = np.exp(-a * T)
    exp_bT  = np.exp(b  * T)
    exp_cT  = np.exp(-c * T)


    C1, C2, C3 = solve_constants(eq_params, ab_params, stage="pre")
    E1, E2, _  = solve_constants(eq_params, ab_params, stage="post")
    
    # consimption rate prop to \dot{g}(t)
    gT0 = C1*exp_bT0 + C2*exp_aT0 + C3*exp_cT0
    gT  = E1  + E2  + C3*exp_cT

    Ts   = (1 / b) * np.log(S0 / C1)
    Ts_p = (1 / b) * np.log((S0 + gT - gT0) / E1) + T

    # consumption rate prop to g(t)
    #G0  = (C1/b) - (C2/a) - (C3/c)
    #GT  = (E1/b) - (E2/a) - (C3/c)*exp_cT
    #GT0 = (C1/b)*exp_bT0 - (C2/a)*exp_aT0 - (C3/c)*exp_cT0

    # consumption time without antibiotics
    #Ts   = (1 / b) * np.log((b / C1)*(S0/(1-δ) + G0))

    # consumption time with antibiotics
    #Ts_p = (1 / b) * np.log((b / E1)*(S0/(1-δ) + G0 - GT0 + GT)) + T

    # solving witout approx.
    # #B,_,_ = solve_constants(eq_params, ab_params, 'post')
    # Ts_min = 0 # np.log(1/f)
    # Ts_max = 100 #np.max(T0 * 40)
    # #Ts_arr = np.linspace(Ts_min, Ts_max, 100)

    # Ts, Ts_p = estimate_Ts(eq_params, ab_params, iter=3)

    # weighted average
    Ts_avrg = (1 - p) * Ts + p * Ts_p

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


def transform_fitness_to_bethedge_parameters(fitness, bac_params):
    '''
    This function identifies the optimal parameters from the maximal fitness
    '''
    x  = bac_params['x']
    λr = bac_params['λr']
    δ  = bac_params['δ']

    F_max = np.max(fitness)
    idx = (fitness == F_max)

    return F_max, x[idx], λr[idx], δ[idx]
