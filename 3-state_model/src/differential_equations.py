# Coupled ODEs that make up the model
import numpy as np

##########################
## Intrinsic parameters ##
##########################
K = 10**9           # carrying capacity
β_max = 1           # max growth rate
γ = 1               # kill rate
λ_min = 0.01        # lower limit on lag time
δ_max = 0.1         # upper limit on δ



#########################
## External parameters ##
#########################
f = 10**(-6)        # dilution fraction
S0 = K              # initial substrate
p0 = f * S0         # initial population
p_min = 1           # lower threshold

dt_max = 1 #check value
T0_max = 12 # check value
Ts_max = 50


############################
## Differential equations ##
############################
# Growth rate as function of nutrients
def β(s):
    return β_max				# constant
    # return beta_max * (s / (s + K))		# Monod


# ODE without antibiotics
def ode_grow_single(t, p, λ_d, λ_r, δ):
    p[p < 0] = 0                                    # avoid negative populations

    dp_dt = np.zeros_like(p)
    dp_dt[0] = -p[0] / λ_d
    dp_dt[1] =  p[0] / λ_d + p[2] / λ_r + p[1] * (1-δ)
    dp_dt[2] = -p[2] / λ_r + p[1] * δ

    return dp_dt


def ode_grow(t, p, λ_d, λ_r, δ):
    p[p < 0] = 0                                    # avoid negative populations

    dp_dt = np.zeros_like(p)
    dp_dt[0:2] = -p[0:2] / λ_d
    dp_dt[2:4] =  p[0:2] / λ_d + p[4:6] / λ_r - p[2:4] * δ
    dp_dt[4:6] = -p[4:6] / λ_r + p[2:4] * δ

    if p[6] > 0:
        dp_dt[2:4] += β(p[6]) * p[2:4]              # adding growth if nutrients left
        #dp_dt[6] = -dp_dt[2:4].sum()               # nutrients, old consumption rate
        dp_dt[6] = -((1-δ) * p[2:4]).sum()          # correct consumption rate


    return dp_dt


# ODE with antibiotics
def ode_kill_single(t, p, λ_d, λ_r, δ):
    p[p < 0] = 0                                    # avoid negative populations

    dp_dt = np.zeros_like(p)
    dp_dt[0] = -p[0] / λ_d
    dp_dt[1] =  p[0] / λ_d + p[2] / λ_r - p[1] * (γ + δ)
    dp_dt[2] = -p[2] / λ_r + p[1] * δ

    return dp_dt


def ode_kill(t, p, λ_d, λ_r, δ):
    p[p < 0] = 0                                    # avoid negative populations

    dp_dt = np.zeros_like(p)
    dp_dt[0:2] = -p[0:2] / λ_d
    dp_dt[2:4] =  p[0:2] / λ_d + p[4:6] / λ_r - p[2:4] * (γ + δ)
    dp_dt[4:6] = -p[4:6] / λ_r + p[2:4] * δ

    return dp_dt
