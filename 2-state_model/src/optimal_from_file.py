import numpy as np
from analytical_calculations import compute_a_and_b, compute_ap_and_bp

# Function that yields optimal parameters of single species winner
def optimal_parameters_from_data(bac_args, ab_args):

    lag, delta, a, b, ap, bp = bac_args
    p, T0, Tab = ab_args
    ab_res = len(np.loadtxt(f"data/low_resolution/optimal_lag-T0{int(T0)}.txt"))
    #ab_res = len(np.loadtxt(f"data/low_resolution/optimal_lag-Tab{int(Tab)}.txt"))
    ip = int(p * (ab_res-1))
    it = int(Tab * (ab_res-1) / 24)

    lag[0]   = np.loadtxt(f"data/low_resolution/optimal_lag-T0{int(T0)}.txt")[ip, it]
    delta[0] = np.loadtxt(f"data/low_resolution/optimal_delta-T0{int(T0)}.txt")[ip, it]
    #lag[0]   = np.loadtxt(f"data/low_resolution/optimal_lag-Tab{int(Tab)}.txt")[ip, it]
    #delta[0] = np.loadtxt(f"data/low_resolution/optimal_delta-Tab{int(Tab)}.txt")[ip, it]

    # Transforming to a-b scheme
    a[0],  b[0]  = compute_a_and_b(lag[0],   delta[0])
    ap[0], bp[0] = compute_ap_and_bp(lag[0], delta[0])

    return [lag, delta, a, b, ap, bp]