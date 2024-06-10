import numpy as np

from differential_equations import S0, f
from analytical_calculations import solve_constants

def total_population_T_p(eq_params, ab_params):
    a  = eq_params['a']
    b  = eq_params['b']
    c  = eq_params['c']
    ap = eq_params['ap']
    bp = eq_params['bp']
    
    Tab = ab_params['Tab']
    T   = ab_params['T']

    d_0 = f*S0

    Bp, Ap, Cp = solve_constants(eq_params, ab_params, 'ab')
    Br = Bp * (ap - a*b) / (a*b)
    Ar = Ap * (bp - a*b) / (a*b)
    Cr = (ap+bp-a*b-c) / (a*b)
    Cr -= c*d_0 / (a*b)

    d_T = d_0*np.exp(-c*T)
    g_T = Bp*np.exp(-bp*Tab) + Ap*np.exp(-ap*Tab) + Cp*np.exp(-c*T)
    r_T = Br*np.exp(-bp*Tab) + Ar*np.exp(-ap*Tab) + Cr*np.exp(-c*T)

    p_T = d_T + g_T + r_T

    return d_T, g_T, r_T # p_T



def total_population_T(eq_params, ab_params):
    a  = eq_params['a']
    b  = eq_params['b']
    c  = eq_params['c']
    T   = ab_params['T']

    d_0 = f*S0

    B, A, C = solve_constants(eq_params, ab_params, 'pre')
    Br = B * (1-b)/b
    Ar = -A * (a+1)/a
    Cr = -C * (a+1)*(1-b) / (c-a*b)
    Cr -= c*d_0 / (a*b)

    d_T = d_0*np.exp(-c*T)
    g_T =  B*np.exp(b*T) +  A*np.exp(-a*T) +  C*np.exp(-c*T)
    r_T = Br*np.exp(b*T) + Ar*np.exp(-a*T) + Cr*np.exp(-c*T)

    p_T = d_T + g_T + r_T

    return d_T, g_T, r_T # p_T



def weighted_population_T(eq_params, ab_params):
    p = ab_params['p']
    d_T, g_T, r_T = total_population_T(eq_params, ab_params)
    d_T, g_T, r_T = total_population_T_p(eq_params, ab_params)

    #p_T_avrg = (1-p) * p_T + p * p_Tp

    return d_T, g_T, r_T #p_T_avrg