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
    nb_trig = 3

    # splitting populations in only triggered persistence, and both
    g_trig = p_t[N:N+Nd]
    g_both = p_t[N+Nd:2*N].reshape(Nδ-1, Nr, Nd)    			# isolating growing populations. Excludin

    # only triggered persistence
    g_pad = np.pad(g_trig, 1, mode='edge')
    λd_mut   = (np.roll(g_pad, 1, axis=0) + np.roll(g_pad, -1, axis=0))[1:Nd+1]
    δ_λr_mut = p_t[N+Nd:N+2*Nd]

    mutation_from_trig = nb_trig * p_t[N:N+Nd]
    mutation_to_trig   = λd_mut + δ_λr_mut

    # both triggered and spontaneous
    # mutation between nearest neighbours in parameters space
    g_pad = np.pad(g_both, 1, mode='edge')

    λd_mut = (np.roll(g_pad, 1, axis=2) + np.roll(g_pad, -1, axis=2))[1:Nδ, 1:Nr+1, 1:Nd+1]
    λr_mut = (np.roll(g_pad, 1, axis=1) + np.roll(g_pad, -1, axis=1))[1:Nδ, 1:Nr+1, 1:Nd+1]
    δ_mut  = (np.roll(g_pad, 1, axis=0) + np.roll(g_pad, -1, axis=0))[1:Nδ, 1:Nr+1, 1:Nd+1]

    mutation_from_both = nb * p_t[N+Nd:2*N]
    mutation_to_both   = (λd_mut + δ_mut + λr_mut).flatten()

    # adding mutation to only triggered persistence
    mutation_from_both[:Nd] = p_t[N+Nd:N+2*Nd]
    mutation_to_both[:Nd]   = p_t[N:N+Nd]

    mutate_from = np.concatenate([mutation_from_trig, mutation_from_both])
    mutate_to   = np.concatenate([mutation_to_trig,   mutation_to_both])

    return mutate_to - mutate_from


# 3 state model of triggered and spontaneous persistence
def ode_grow(t, p_t, λd, λr, δ, Ɛ):
    p_t[:-1][p_t[:-1] < 0] = 0 #10**(-200)     			# avoid negative populations

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
        dn_dt[N:2*N]   =  p_t[:N] / λd + p_t[2*N:3*N] / λr  + p_t[N:2*N] * (1 - δ) # - nb*Ɛ)
        dn_dt[2*N:3*N] = -p_t[2*N:3*N] / λr + p_t[N:2*N]*δ

        # mutation
        dn_dt[N:2*N] += Ɛ * mutate(p_t, λd, λr, δ)

        # update nutrients
        dn_dt[-1] = -dn_dt[N:2*N].sum()

    return dn_dt


def ode_kill(t, p_t, λd, λr, δ, Ɛ):
    # avoid negative populations
    p_t[:-1][p_t[:-1] < 0] = 0 # 10 ** (-200)

    # number of species
    N = len(λd)
    
    dn_dt = np.zeros_like(p_t)
    dn_dt[0:N]     = -p_t[:N] / λd
    dn_dt[N:2*N]   =  p_t[:N] / λd + p_t[2*N:3*N] / λr - p_t[N:2*N] * (γ + δ)
    dn_dt[2*N:3*N] = -p_t[2*N:3*N] / λr + p_t[N:2*N]*δ 

    return dn_dt

