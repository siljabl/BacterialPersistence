import numpy as np

from differential_equations import S0
from analytical_calculations import solve_constants


def compute_substrate(t, eq_params, ab_params):
    '''
    Computing amount of substrate at time t, for cycles without antibiotics.
    The amount of substrate is equal to the initial amount (S0) minus the size of the growing population (g_t).
    t: time
    '''
    a  = eq_params['a']
    b  = eq_params['b']
    c  = eq_params['c']

    B0, A0, C0 = solve_constants(eq_params, ab_params, 'pre')

    g_t = B0*np.exp(b*t) + A0*np.exp(-a*t) + C0*np.exp(-c*t)

    return S0 - g_t


def compute_substrate_p(t, eq_params, ab_params):
    '''
    Computing amount of substrate at time t, for cycles with antibiotics.
    The amount of substrate is equal to the initial amount (S0) minus the size of the growing population (g_t - g_T + g_T0).
    t: time
    '''
    a  = eq_params['a']
    b  = eq_params['b']
    c  = eq_params['c']
    ap = eq_params['ap']
    bp = eq_params['bp']
    
    T0  = ab_params['T0']
    Tab = ab_params['Tab']
    T   = ab_params['T']

    B0, A0, C0 = solve_constants(eq_params, ab_params, 'pre')
    Bp, Ap, Cp = solve_constants(eq_params, ab_params, 'ab')
    B,  A,  C  = solve_constants(eq_params, ab_params, 'post')
    
    g_T0 = B0*np.exp( b*T0)   + A0*np.exp(-a*T0)   + C0*np.exp(-c*T0)
    g_T  = Bp*np.exp(-bp*Tab) + Ap*np.exp(-ap*Tab) + Cp*np.exp(-c*T)
    g_t  =  B*np.exp( b*t)    +  A*np.exp(-a*t)    +  C*np.exp(-c*t)

    return S0 - (g_t - g_T + g_T0)
