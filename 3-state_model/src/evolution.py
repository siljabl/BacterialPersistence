import numpy as np
from tqdm import tqdm
from scipy.integrate import solve_ivp

from differential_equations import S0, f, p_min, dt_max, Ts_max
from mutation_differential_equations import ode_grow, ode_kill


#######################
## Solving one cycle ##
#######################
def total_population(p_t, N):
    '''
    Computing the total population of each species, by summing over the dormant (d_t), persistent (r_t) and growing (g_t) subpopulations
    '''
    d_t = p_t[:N]
    g_t = p_t[N:2*N]
    r_t = p_t[2*N:3*N]

    p_sum = d_t + g_t + r_t

    return p_sum



def return_to_sleep(p_t, N):
    '''
    Putting the total population of all species to sleep
    '''
    d_t = total_population(p_t, N)
    g_t = np.zeros(N)
    r_t = np.zeros(N)
    p_t = np.concatenate([d_t, g_t, r_t, np.array([S0])])

    return p_t



def remove_extinct_species(p_t, mask_alive, N):
    '''
    Removing an extinct species by setting all subpopulations to 0.
    A species is considered extinct it its total population is below the threshold p_min.
    '''
    p_sum = total_population(p_t, N)
            
    p_t[0*N:1*N][p_sum < p_min] = 0
    p_t[1*N:2*N][p_sum < p_min] = 0
    p_t[2*N:3*N][p_sum < p_min] = 0
        
    n_dead = sum((p_sum < p_min) * mask_alive)	# updating counter
    mask_alive = mask_alive & (p_sum > p_min)	# updating mask to contain only species that are alive after ab

    return n_dead, mask_alive



def famine(t, p_t, lag, delta, omega, eps): return p_t[-1]
'''
Function for detecting when system runs out of food. Used by solve_ivp().
'''
famine.terminal  = True
famine.direction = -1


def sort_sol_cycle(sol_cycle):
    if len(sol_cycle) == 3:
        substrate = np.concatenate([sol_cycle[0].y[-1], sol_cycle[1].y[-1], sol_cycle[2].y[-1]])
        species = np.concatenate([sol_cycle[0].y[:-1], sol_cycle[1].y[:-1], sol_cycle[2].y[:-1]], axis=1)
        time = np.concatenate([sol_cycle[0].t, sol_cycle[1].t, sol_cycle[2].t])

    else:
        substrate = sol_cycle.y[-1]
        species = sol_cycle.y[:-1]
        time = sol_cycle.t

    return [time, species, substrate]



def solve_cycle(p_t, r, bac_params, ab_params, sim_params):
    '''
    Performing one feast-famine cycle. The final (sub)populations are found with solve_ivp.
    The number of species that have gone extinct are counted after every antibiotic or dillusion event.
    '''
    p  = ab_params['p']
    T0 = ab_params['T0']
    T  = ab_params['T']

    check_extinction = sim_params['extinct']
    N = len(bac_params['λd'])
    n_states = 3

    # creating mask of none-extinct species before antibiotic is applied
    p_sum      = total_population(p_t, N)
    mask_alive = (p_sum >= p_min)

    Ts_min  = 0
    n_dead = 0

    args = [bac_params['λd'], bac_params['λr'], bac_params['δ'], bac_params['Ɛ']]

    # cycles with antibiotics
    if r < p:
        # solving for period before antibiotics
        sol_cycle_1 = solve_ivp(ode_grow, [0, T0], p_t, args=args, max_step=dt_max)
        p_t = np.array([sol_cycle_1.y[i][-1] for i in range(n_states * N + 1)])

        # solving for period during antibiotics
        sol_cycle_2 = solve_ivp(ode_kill, [T0, T], p_t, args=args, max_step=dt_max)
        p_t = np.array([sol_cycle_2.y[i][-1] for i in range(n_states * N + 1)])
        
        # checking for n_extinct (only for simulations with lower threshold implimented)
        if check_extinction == True:
            n_dead, mask_alive = remove_extinct_species(p_t, N)

        # updating lower bound for integration
        Ts_min = T

    # cycles without antibiotics / solving for period after antibiotics
    sol_cycle = solve_ivp(ode_grow, [Ts_min, Ts_max], p_t, args=args, events=famine, max_step=dt_max)
    p_t = np.array([sol_cycle.y[i][-1] for i in range(n_states * N + 1)])
   
    p_t = return_to_sleep(p_t, N)
    p_t[:N] = f * p_t[:N]

    # cheking for n_extinct (only for simulations with lower threshold implimented)
    if check_extinction == True:
        n_dead += remove_extinct_species(p_t, mask_alive, N)[0]

    # sort output of solve_ivp for plotting
    if r < p:
        sol_cycle = [sol_cycle_1, sol_cycle_2, sol_cycle]
    [time, species, substrate] = sort_sol_cycle(sol_cycle)
        
    return p_t, n_dead, [time, species, substrate]



def evolve_system(p_t, bac_params, ab_params, sim_params):
    '''
    Solving several feast-famine cycles to find optimal parameters after N cycles.
    '''
    # breaking if mutation rate = 0
    λd = bac_params['λd']
    λr = bac_params['λr']
    δ  = bac_params['δ']
    N  = len(λd)

    # sorting parameters
    r_arr      = sim_params['r_arr']
    tot_cycles = sim_params['tot_cycles']

    # output arrays
    λd_avrg = np.zeros_like(r_arr)
    λr_avrg = np.zeros_like(r_arr)
    δ_avrg  = np.zeros_like(r_arr)

    p_dominant = np.zeros_like(r_arr)
    n_extinct  = np.zeros_like(r_arr)
    p_dists    = np.zeros([len(r_arr), N])

    species = []
    substrate = []
    time = []

    # looping through cycles
    for ic in tqdm(range(tot_cycles)):
        p_t, n_dead, sol_cycle = solve_cycle(p_t, r_arr[ic], bac_params, ab_params, sim_params)
        substrate.append(sol_cycle[2])
        species.append(sol_cycle[1][:N] + sol_cycle[1][N:2*N] + sol_cycle[1][2*N:3*N])
        time.append(sol_cycle[0])
        
        d_t = p_t[:N]
        # if all species go extinct
        if sum(d_t) == 0:
            break

        λd_avrg[ic] = sum(λd * d_t) / sum(d_t)
        λr_avrg[ic] = sum(λr * d_t) / sum(d_t)
        δ_avrg[ic]  = sum(δ  * d_t) / sum(d_t)

        p_dominant[ic] = np.where(d_t == d_t.max())[0][0]
        n_extinct[ic]  = n_dead
        p_dists[ic]    = d_t / sum(d_t)
        
    return λd_avrg, λr_avrg, δ_avrg, p_dominant, n_extinct, p_dists, [species, substrate, time]
   
