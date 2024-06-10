import numpy as np

from differential_equations import γ


############################
## Differential equations ##
############################
def mutate(p_t, λd, λr, δ):
    Nd = len(np.unique(λd))
    Nr = len(np.unique(λr))
    Nδ = len(np.unique(δ))
    N  = len(λd)
    nb = 6
        
    g_t = p_t[N:2*N].reshape(Nδ, Nd, Nr)    			# isolating growing populations

    # mutation between nearest neighbours in parameters space
    g_pad = np.pad(g_t, 1, mode='edge')

    λd_mut = (np.roll(g_pad, 1, axis=1) + np.roll(g_pad, -1, axis=1))[1:Nr+1, 1:Nd+1, 1:Nδ+1]
    λr_mut = (np.roll(g_pad, 1, axis=0) + np.roll(g_pad, -1, axis=0))[1:Nr+1, 1:Nd+1, 1:Nδ+1]
    δ_mut  = (np.roll(g_pad, 1, axis=2) + np.roll(g_pad, -1, axis=2))[1:Nr+1, 1:Nd+1, 1:Nδ+1]

    #mutate_away = nb * p_t[N:2*N]
    mutate_to = (λd_mut + δ_mut + λr_mut).flatten()

    return mutate_to #- mutate_away


# 3 state model of triggered and spontaneous persistence
def ode_grow(t, p_t, λd, λr, δ, Ɛ):
    p_t[:-1][p_t[:-1] < 0] = 10**(-200)     			# avoid negative populations

    Nd = len(np.unique(λd))
    Nr = len(np.unique(λr))
    Nδ = len(np.unique(δ))
    N = len(λd)
    nb = 6

    dn_dt = np.zeros_like(p_t)              			# empty array for computing derivatives
    # g_t = p_t[N:2*N].reshape(Nδ, Nd, Nr)
    # all mutate normally, except species with only triggered persistence
    # they mutate in both wake up rate and persistence rate

    if p_t[-1] > 0:
        # normal differential equations
        dn_dt[:N]      = -p_t[:N] / λd
        dn_dt[N:2*N]   =  p_t[:N] / λd + p_t[2*N:3*N] / λr  + p_t[N:2*N] * (1 - δ)# - nb*Ɛ)
        dn_dt[2*N:3*N] = -p_t[2*N:3*N] / λr + p_t[N:2*N]*δ

        # mutation
        #dn_dt[N:2*N] += Ɛ * mutate(p_t, λd, λr, δ)

        # update nutrients
        dn_dt[-1] = -dn_dt[N:2*N].sum()

    return dn_dt


def ode_kill(t, p_t, λd, λr, δ, Ɛ):
    # avoid negative populations
    p_t[:-1][p_t[:-1] < 0] = 10 ** (-200)

    # number of species
    N = len(λd)
    
    dn_dt = np.zeros_like(p_t)
    dn_dt[0:N]     = -p_t[:N] / λd
    dn_dt[N:2*N]   =  p_t[:N] / λd + p_t[2*N:3*N] / λr - p_t[N:2*N] * (γ + δ)
    dn_dt[2*N:3*N] = -p_t[2*N:3*N] / λr + p_t[N:2*N]*δ 

    return dn_dt

