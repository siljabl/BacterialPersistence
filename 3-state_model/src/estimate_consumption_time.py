import numpy as np

from differential_equations import S0, f
from substrate import compute_substrate, compute_substrate_p


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
    args['idx_sort'] = np.lexsort((args['idx'][3], args['idx'][2], args['idx'][1]))

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
    S2 = S[args['idx'][0] + 1, args['idx'][1], args['idx'][2], args['idx'][3]]
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


def define_initial_Ts_arr(ab_params):
    Ts_min = 10
    Ts_max = 24
    Ts_p_max = 4 * ab_params['Tab']
    if Ts_p_max < Ts_max:
        Ts_p_max = Ts_max
    
    Ts_arr = np.linspace(Ts_min,   Ts_max, 100)
    Ts_p_arr = np.linspace(Ts_min, Ts_p_max, 100)

    return Ts_arr, Ts_p_arr



def estimate_Ts(eq_params, ab_params, iter=5):
    '''
    This function estimates the time at which all substrate is consumed (Ts).

    Ts_arr: range of Ts that is explored.
    iter: number of iterations for improving precision.
    '''
    dims  = np.shape(eq_params['a'])
    NA    = np.newaxis

    Ts_arr, Ts_p_arr = define_initial_Ts_arr(ab_params)
    steps = len(Ts_arr)

    args   = {'Ts':Ts_arr,   'dt2':Ts_arr[1]   - Ts_arr[0]}
    args_p = {'Ts':Ts_p_arr, 'dt2':Ts_p_arr[1] - Ts_p_arr[0]}

    args['Ts_min']   = (np.random.rand(*dims) - 0.5) * args['dt2']
    args_p['Ts_min'] = (np.random.rand(*dims) - 0.5) * args['dt2']


    # increasing precision on Ts
    for i in range(iter):
        args['temp']   = args['Ts'][:, NA, NA, NA]   + args['Ts_min']
        args_p['temp'] = args_p['Ts'][:, NA, NA, NA] + args_p['Ts_min']
        
        # analytical population
        S  = compute_substrate(args['temp'],   eq_params, ab_params)
        Sp = compute_substrate_p(args_p['temp'], eq_params, ab_params)

        # Finding upper and lower limit on ts
        update_Ts_min(S,  args,   dims)
        update_Ts_min(Sp, args_p, dims)

        # Preparing for next round
        update_time_arrays(args, steps)
        update_time_arrays(args_p, steps)


    # Linear approximation around S = 0
    S1, S2 = reshape_S_arrays(S, args, dims)
    Ts = linear_interpolation(args['dt2'], args['Ts'], S1, S2)

    S1, S2 = reshape_S_arrays(Sp, args_p, dims)
    Ts_p = linear_interpolation(args_p['dt2'], args_p['Ts'], S1, S2)

    return Ts, Ts_p


