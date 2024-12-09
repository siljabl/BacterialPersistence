# Functions for computing optimal parameters from functions in model_equations.py
import numpy as np

from differential_equations  import S0, f, n_min
from analytical_calculations import analytical_growth, analytical_decay


#################
## Competition ##
#################

def update_Ts_min(S, args, dims):
    '''
    This function identifies the indices in the Ts-array before and after the substrate is consumed (i.e. before and after S intersects with 0).
    These are then saved as 'Ts_min' to be used as new boundaries on Ts.
    
    S: substrate, 3-dim array
    args: arguments used to estimate Ts
    dims: dimensions of S-array
    '''

    # find indices where S changes sign and rearrange indices to 3-dim array (like S)
    S_diff = np.diff(np.sign(S), axis=0)
    args['idx'] = np.asarray(S_diff < 0).nonzero()
    args['idx_sort'] = np.lexsort((args['idx'][2], args['idx'][1]))

    # use index of S just before crossing 0 as lower limit on Ts
    lowlim = args['Ts'][args['idx'][0]]
    lowlim = lowlim[args['idx_sort']].reshape(dims) + args['Ts_min']
    args['Ts_min'] = lowlim


def update_time_arrays(args, t_len):
    '''
    This function updates the values of Ts, dt1, dt2, for a new iteration of the main function.

    args: arguments used to estimate Ts
    t_len: length of Ts array
    '''
    # update Ts array with improved resolution
    args['Ts']  = np.linspace(0, args['dt2'], t_len)
    
    # saving old resolution
    args['dt1'] = args['dt2']

    # updating new resolution
    args['dt2'] = args['Ts'][1] - args['Ts'][0]



def reshape_S_arrays(S, args, dims):
    '''
    This function defines S1 and S2, i.e. the value of the substrate before and after it becomes zero, and reshapes the array to the shape of S.
    
    '''
    # saving t just before intersection with 0 as Ts, to use in linear interpolation
    args['Ts'] = args['Ts_min']
    
    # saving S just before intersection as S1
    S1 = S[args['idx']]
    S1 = S1[args['idx_sort']].reshape(dims)

    # saving S just after intersection as S2
    S2 = S[args['idx'][0] + 1, args['idx'][1], args['idx'][2]]
    S2 = S2[args['idx_sort']].reshape(dims)

    return S1, S2


def linear_interpolation(dx, x1, y1, y2):
    '''
    Interpolate to find x: y(x)=0
    '''
    dy = y2 - y1
    a = dy / dx
    b = y1 - a*x1

    return -b / a


# Function for finding ts = t: S(t) = 0
def Ts_approximation(t, n_t0, a, b, iter=3):
    dims  = np.shape(a[0])
    d_t0, g_t0, S_0 = n_t0                                          # initial values
    form = np.shape(S_0)                                            # shape of output
    NA = np.newaxis                                                 # extra axis for array multiplication

    t_steps = len(t)                                                # length of time array

    args = {'Ts':t, 'dt2':t[1]-t[0]}
    args['Ts_min'] = (np.random.rand(form[0], form[1]) - 0.5) * args['dt2'] 

    # increasing precision on Ts
    for i in range(iter):
        args['temp'] = args['Ts'][:, NA, NA, NA] + args['Ts_min']

        #t_temp = t[:, NA, NA, NA] + t_rand                          # temporary time array
        exp_bt = np.exp(b * args['temp'])                                 # time dependent terms
        exp_at = np.exp(-a * args['temp'])                                # time dependent terms
        
        # analytical population
        g_t = d_t0 * a * b * (exp_bt - exp_at) / (a + b) + \
              g_t0 * (b * (1 + a) * exp_bt + a * (1 - b) * exp_at) / (a + b)  	# size of growing population
        S = S_0 - (g_t - g_t0).sum(axis=len(np.shape(g_t)) - 2 - 1)             # amount of substrate
        #print("S:", np.shape(S))


        update_Ts_min(S, args, dims)
        update_time_arrays(args, t_steps)
   
    S1, S2 = reshape_S_arrays(S, args, dims)
    Ts = linear_interpolation(args['dt2'], args['Ts'], S1, S2)

    return Ts



# Finding optimal competitor parameters
def run_competition(bac_args, ab_args, sim_args):
    # Sorting inputs
    lag, delta, a, b, ap, bp = bac_args
    p, T0, Tab = ab_args[0:3]
    _, bac_res, t_res, tot_cycles = sim_args

    t_min, t_max = 10 - T0, 15 + (T0 + Tab)                     # time limits
    t_arr = np.linspace(t_min, t_max, t_res)                    # time array

    S_frac_cycle = np.zeros([tot_cycles, 2, bac_res, bac_res])  # array for output


    # Initial populations, as [d(0), g(0), S(0)]
    n_0 = [f * S0 * np.ones_like(lag), np.zeros_like(lag), S0 * np.ones_like(lag)[0]]
    extinct = np.zeros([bac_res, bac_res])                    # array for counting first extinctions
    ext     = np.ones_like(lag)                               # masking 
    r_arr   = np.random.rand(tot_cycles)                      # random array for ab

    for ic in range(tot_cycles):

        # With antibiotics
        if r_arr[ic] < p:
            n_T0 = analytical_growth(T0, n_0,  a, b)            # population before AB
            n_T  = analytical_decay(Tab, n_T0, a, b, ap, bp)    # population after AB

            d_dead = n_T[0] < n_min                             # checking if dormant species is killed by AB
            g_dead = n_T[1] < n_min                             # checking if growing species is killed by AB
            ext[d_dead * g_dead * (extinct == 0)] = 0
            for i in range(2):
                extinct[d_dead[i] * g_dead[i]] = 1              # updating counter of first extinction

            t   = Ts_approximation(t_arr, n_T, a, b)            # computing ts: S(ts) = 0
            n_t = analytical_growth(t, n_T, a, b)               # population at ts
            S_frac_cycle[ic] += (n_t[1] - n_T[1] + n_T0[1]) / (S0 - n_t[2])

        # Without antibiotics
        else:
            t   = Ts_approximation(t_arr, n_0, a, b)            # computing ts: S(ts) = 0
            n_t = analytical_growth(t, n_0, a, b)               # population at ts
            S_frac_cycle[ic] += n_t[1] / (S0 - n_t[2])

        # Preparing for next cycle
        d_0  = f * (n_t[0] + n_t[1])                            # enter dormancy and dilute
        dead = d_0 < n_min                                      # checking if species killed by dilution
        ext[dead * (extinct == 0)] = 0
        for i in range(2):
            extinct[dead[i]] = 1                                # updating counter of first extinction
        n_0 = [d_0 * np.ones_like(lag), np.zeros_like(lag), S0 * np.ones_like(lag)[0]]

    # Finding optimal set of bacterial parameters
    S_frac_mean = S_frac_cycle[:, 1].mean(axis=0)                   # taking cycle average of consumption fraction
    S_max = (S_frac_mean == S_frac_mean.max())                      # finding max consumption fraction

    S_frac  = S_frac_mean.max()                                     # saving max consumption fraction
    lag_opt = lag[1,0][(S_max.sum(0)).astype(bool)][0]              # saving corresponding lag time
    del_opt = delta[1][(S_max.sum(1)).astype(bool), 0][0]           # saving corresponding type-II fraction


    return np.array([S_frac, lag_opt, del_opt]), n_t, t, S_frac_cycle, r_arr
