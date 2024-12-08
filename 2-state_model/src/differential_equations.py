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
f = 2*10**(-6)        # dilution fraction
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


