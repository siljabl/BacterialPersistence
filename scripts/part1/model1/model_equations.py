# Coupled ODEs that make up the model
import numpy as np

##########################
## Intrinsic parameters ##
##########################
K = 10**9				# carrying capacity
beta_max = 1			# max growth rate
gamma = 1 				# kill rate
lag_min = 0.01			# lower limit on lag time


#########################
## External parameters ##
#########################
f = 10**(-6)       	# dilution fraction
S0 = K             	# initial substrate
n0 = f * S0         	# initial population
n_min = 1000        	# lower threshold
T_cycle = 50		# upper limit on cycle length



############################
## Differential equations ##
############################
# Growth rate as function of nutrients
def beta(s):
    return beta_max				# constant
    # return beta_max * (s / (s + K))		# Monod


# Computing derivatives in absence of antibiotics
def growth(t, n, lag):
    n[:-1][n[:-1] < 0] = 10**(-200)     		# avoid negative populations

    # array of derivatives
    dn_dt = np.zeros_like(n)
    dn_dt[0:2] = -n[0:2] / lag			# dormant populations
    dn_dt[2:4] = n[0:2] / lag			# growing populations

    # growth phase
    if n[4] > 0:
        dn_dt[2:4] += beta(n[4]) * n[2:4]	# growing populations
        dn_dt[4] = -sum(dn_dt[2:4])		# nutrients

    return dn_dt


# Computing derivatives in presence of antibiotics
def decay(t, n, lag):
    n[:-1][n[:-1] < 0] = 10**(-200)     		# avoid negative populations

    # array of derivatives
    dn_dt = np.zeros_like(n)
    dn_dt[0:2] = -n[0:2] / lag			# dormant populations
    dn_dt[2:4] = n[0:2] / lag - gamma * n[2:4]	# growing populations

    return dn_dt

