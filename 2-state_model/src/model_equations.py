# Coupled ODEs that make up the model
import numpy as np

##########################
## Intrinsic parameters ##
##########################
K = 10**9           # carrying capacity
beta_max = 1        # max growth rate
gamma = 1           # kill rate
lag_min = 0.01      # lower limit on lag time
delta_max = 0.1     # upper limit on delta


#########################
## External parameters ##
#########################
f = 10**(-6)        # dilution fraction
S0 = K              # initial substrate
n0 = f * S0         # initial population
n_min = 1           # lower threshold


############################
## Differential equations ##
############################
# Growth rate as function of nutrients
def beta(s):
    return beta_max				# constant
    # return beta_max * (s / (s + K))		# Monod


# ODE without antibiotics
def ode_grow(t, n, lag, delta):
    n[n < 0] = 0                                    # avoid negative populations

    dn_dt = np.zeros_like(n)
    dn_dt[0:2] = -n[0:2] / lag + delta * n[2:4]     # dormant populations
    dn_dt[2:4] = -dn_dt[0:2]                        # awake populations

    if n[4] > 0:
        dn_dt[2:4] += beta(n[4]) * n[2:4]           # adding growth if nutrients left
        dn_dt[4] = -dn_dt[2:4].sum()                # nutrients

    return dn_dt


# ODE with antibiotics
def ode_kill(t, n, lag, delta):
    n[n < 0] = 0                                    # avoid negative populations

    dn_dt = np.zeros_like(n)
    dn_dt[0:2] = -n[0:2] / lag + delta * n[2:4]     # dormant populations
    dn_dt[2:4] = -dn_dt[0:2] - gamma * n[2:4]       # awake populations

    return dn_dt


##########################
## a-b scheme equations ##
##########################
# computing a, b, ap, bp as defined in thesis
def a_b(lag, delta):
    D = np.sqrt((lag * (1 - delta) - 1) ** 2 + 4 * lag)

    a = (D - (lag * (1 - delta) - 1)) / (2 * lag)
    b = (D + (lag * (1 - delta) - 1)) / (2 * lag)

    return a, b


def ap_bp(lag, delta):
    D = np.sqrt((lag * (gamma + delta) + 1) ** 2 - 4 * gamma * lag)

    ap = (D + (lag * (gamma + delta) + 1)) / (2 * lag)
    bp = - (D - (lag * (gamma + delta) + 1)) / (2 * lag)

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
