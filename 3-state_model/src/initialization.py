import numpy as np

from differential_equations import S0, f, λ_min, δ_max

###########################
## Simulation parameters ##
###########################
def initialize_bacterial_parameter_arrays(ab_params, dλ, dδ):
    T   = ab_params['T']

    # determining size of system
    Nd = int(T / dλ) + 1
    Nr = int(T / dλ) + 1
    Nδ = int(δ_max / dδ) + 1

    # defining parameter arrays
    λd = dλ * np.arange(0, Nd, 1) + λ_min
    λr = dλ * np.arange(0, Nr, 1) + λ_min
    δ  = dδ * np.arange(0, Nδ, 1)

    λd = np.tile(np.tile(λd, Nδ), Nr)
    λr = np.tile(np.repeat(λr, Nd), Nδ)
    δ  = np.repeat(np.repeat(δ, Nd), Nr)
    
    mask_identical_species = np.logical_not((δ == 0) * (λr > λ_min))

    λd = λd[mask_identical_species]
    λr = λr[mask_identical_species]
    δ  = δ[mask_identical_species]

    bac_params = {'λd':λd, 'λr':λr, 'δ':δ, 'Nd':Nd, 'Nr':Nr, 'Nδ':Nδ}

    return bac_params


def optimal_parameters(bac_params, ab_params, path):
    p   = ab_params['p']
    T0  = ab_params['T0']
    T   = ab_params['T']
    T0_max = ab_params['T0_max']

    # determining size of system
    Nd = bac_params['Nd']
    Nr = bac_params['Nr']
    Nδ = bac_params['Nδ']
    N = Nd * Nδ * Nr

    # defining parameter arrays
    λd = bac_params['λd']
    λr = bac_params['λr']
    δ  = bac_params['δ']

    ab_res = len(np.loadtxt(f'../data/optimal_λd-Tab_{int(T)}'))

    ip = int(p * ab_res)
    it = int(T0 * ab_res / T0_max)
    λd_opt = np.loadtxt(f'../data/optimal_λd-Tab_{int(T)}')[ip, it]
    λr_opt = np.loadtxt(f'../data/optimal_λr-Tab_{int(T)}')[ip, it]
    δ_opt  = np.loadtxt(f'../data/optimal_δ-Tab_{int(T)}')[ip, it]

    # corresponding position in parameter space
    iλd = np.where(abs(λd - λd_opt) == min(abs(λd - λd_opt)))[0][0]
    iλr = np.where(abs(λr - λr_opt) == min(abs(λr - λr_opt)))[0][0]
    iδ  = np.where(abs(δ - δ_opt) == min(abs(δ - δ_opt)))[0][0]

    #index = iδ + iλd + il * Nd * Nδ
    index = iδ + iλd + iλr

    return index



def initialise_system(bac_params, sim_params):
    mutation = sim_params['mutation']
    mut_seed = sim_params['mut_seed']

    # determining size of system
    λd = bac_params['λd']
    Nd = bac_params['Nd']
    N  = len(λd)

    d0  = f*S0
    p_t = np.concatenate([np.zeros(3*N), np.array([S0])])

    # choosing whether to allow mutation and initialising initial conditions in parameter space
    if mutation == False: 
        p_t[:N] = d0 / N

    else:        
        # initialising simulation from combination of smallest parameters
        if   mut_seed == 'min':  index = 0
        elif mut_seed == 'max':  index = Nd - 1
        elif mut_seed == 'best': index = optimal_parameters()
        
        p_t[index] = d0
                
    return p_t
 