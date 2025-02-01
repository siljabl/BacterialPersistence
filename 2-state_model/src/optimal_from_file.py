import numpy as np
from analytical_calculations import compute_a_and_b, compute_ap_and_bp

# Function that yields optimal parameters of single species winner
def optimal_parameters_from_data(bac_args, ab_args, const):
    lag, delta, a, b, ap, bp = bac_args
    p, T0, Tab = ab_args

    if const=='Tab':
        lag_data = np.loadtxt(f"data/low_resolution/optimal_lag-Tab{int(Tab)}.txt")
        del_data = np.loadtxt(f"data/low_resolution/optimal_delta-Tab{int(Tab)}.txt")
        ab_res = len(lag_data)
        ip = int(p * (ab_res-1))
        it = int(T0 * (ab_res-1) / 12)

    elif const=='T0':
        lag_data = np.loadtxt(f"data/low_resolution/optimal_lag-T0{int(T0)}.txt")
        del_data = np.loadtxt(f"data/low_resolution/optimal_delta-T0{int(T0)}.txt")
        ab_res = len(lag_data)
        ip = int(p * (ab_res-1))
        it = int(Tab * (ab_res-1) / 24)

    lag[0]   = lag_data[ip, it]
    delta[0] = del_data[ip, it]

    # Transforming to a-b scheme
    a[0],  b[0]  = compute_a_and_b(lag[0],   delta[0])
    ap[0], bp[0] = compute_ap_and_bp(lag[0], delta[0])

    return [lag, delta, a, b, ap, bp]
