import numpy as np
from differential_equations import gamma, n0, f, S0


##########################
## a-b scheme equations ##
##########################
def compute_a_and_b(lag, delta):
    '''
    computing a and b as defined in thesis
    '''
    D = np.sqrt((lag * (1 - delta) - 1) ** 2 + 4 * lag)

    a = (D - (lag * (1 - delta) - 1)) / (2 * lag)
    b = (D + (lag * (1 - delta) - 1)) / (2 * lag)

    return a, b


def compute_ap_and_bp(lag, delta):
    '''
    computing a_p and b_p as defined in thesis
    '''
    D = np.sqrt((lag * (gamma + delta) + 1) ** 2 - 4 * gamma * lag)

    ap =  (D + (lag * (gamma + delta) + 1)) / (2 * lag)
    bp = -(D - (lag * (gamma + delta) + 1)) / (2 * lag)

    return ap, bp


##########################
## Analytical equations ##
##########################
# Single species total population without limitation on food
def analytical_population(t, bac_args, ab_args):
    _, _, a, b, ap, bp = bac_args
    _, T0, Tab = ab_args
    T = T0 + Tab

    # Precomputations for readability
    prefactor = n0 / ((a + b) ** 2 * (ap - bp))
    a_bp = a - bp
    a_ap = a - ap
    b_ap = b + ap
    b_bp = b + bp

    exp_aT0 = np.exp(-a * T0)
    exp_bT0 = np.exp(b * T0)
    exp_apT = np.exp(-ap * (T - T0))
    exp_bpT = np.exp(-bp * (T - T0))

    b_term = (b_ap * (b_bp * exp_aT0 + a_bp * exp_bT0) * exp_bpT -
              b_bp * (b_ap * exp_aT0 + a_ap * exp_bT0) * exp_apT) * np.exp(b * (t - T))
    a_term = (a_bp * (b_ap * exp_aT0 + a_ap * exp_bT0) * exp_apT -
              a_ap * (b_bp * exp_aT0 + a_bp * exp_bT0) * exp_bpT) * np.exp(-a * (t - T))

    return prefactor * (a * b_term + b * a_term)


# Analytical population during growth
def analytical_growth(t, n_t0, a, b):
    d_t0, g_t0, S_t0 = n_t0

    # Precomputations
    denom = a + b
    exp_bt = np.exp(b * t)
    exp_at = np.exp(-a * t)

    # Dormant
    d_t = (d_t0 + g_t0) * (a * (1 - b) * exp_bt + b * (1 + a) * exp_at) / denom + \
          g_t0 * ((1 - b) * exp_bt - (1 + a) * exp_at) / denom
    # Growing
    g_t = a * b * (d_t0 + g_t0) * (exp_bt - exp_at) / denom + \
          g_t0 * (b * exp_bt + a * exp_at) / denom
    # Nutrient, summin along right axis depending on shapes of input
    if len(np.shape(t)) == 1:
        S = S_t0 - (g_t - g_t0).sum(axis=1)
    else:
        S = S_t0 - (g_t - g_t0).sum(axis=len(np.shape(g_t))-2-1)

    return d_t, g_t, S


# Analytical population during antibiotics
def analytical_decay(t, n_t0, a, b, ap, bp):
    d_t0, g_t0, S_t0 = n_t0

    # Precomputations
    denom = ap - bp
    exp_bpt = np.exp(-bp * t)
    exp_apt = np.exp(-ap * t)

    # Dormant
    d_t = d_t0 * ((ap - a * b) * exp_bpt + (a * b - bp) * exp_apt) / denom + \
          g_t0 * (a + 1) * (1 - b) * (exp_bpt - exp_apt) / denom
    # Growing
    g_t = d_t0 * a * b * (exp_bpt - exp_apt) / denom + \
          g_t0 * ((a * b - bp) * exp_bpt - (a * b - ap) * exp_apt) / denom

    return d_t, g_t, S_t0 * np.ones_like(t)



##########################
## Analytical equations ##
##########################
def solve_constants(eq_params, ab_params, stage):
    '''
    H3omputing constants from equations for growing populations.
    '''
    a  = eq_params['a']
    b  = eq_params['b']
    ap = eq_params['ap']
    bp = eq_params['bp']

    T0  = ab_params['T0']
    T   = ab_params['T']

    ebT0, eaT0 = np.exp(b*T0), np.exp(-a*T0)
    ebT,  eaT  = np.exp(-bp*(T-T0)), np.exp(-ap*(T-T0))

    C1 =  f*S0 * (a*b) / (a+b)
    C2 = -f*S0 * (a*b) / (a+b)

    if stage == 'pre':
        return C1, C2
    
    D1 = C1*((a-bp)*ebT0 + (b+bp)*eaT0) / (ap-bp)
    D2 = C2*((a-ap)*ebT0 + (b+ap)*eaT0) / (ap-bp)

    if stage == 'ab':
        return D1, D2
    
    E1 = ((b+ap)*D1*ebT + (b+bp)*D2*eaT) / (a+b)
    E2 = ((a-ap)*D1*ebT + (a-bp)*D2*eaT) / (a+b)

    if stage == 'post':
        return E1, E2,

    else:
        return 0, 0


