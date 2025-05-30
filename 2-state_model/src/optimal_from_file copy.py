import numpy as np


def identify_optimal_parameters_const_T0(ab_params, T0, folder):
    p   = ab_params['p']
    T0  = ab_params['T0']
    Tab = ab_params['Tab']

    λ_file = np.loadtxt(f"{folder}/optimal_lag-T0{T0}.txt")
    δ_file = np.loadtxt(f"{folder}/optimal_delta-T0{T0}.txt")

    # transform p, Tab to indices
    ab_res  = int(400)
    Tab_max = int(24)
    p_arr   = np.linspace(0, 1,       ab_res)
    T_arr   = np.linspace(0, Tab_max, ab_res)

    ip = np.where(abs(p_arr-p)   == np.min(abs(p_arr-p)))
    iT = np.where(abs(T_arr-Tab) == np.min(abs(T_arr-Tab)))

    # find optimal
    λ_opt = λ_file[ip,iT][0][0]
    δ_opt = δ_file[ip,iT][0][0]

    return λ_opt, δ_opt
    

def identify_optimal_parameters_const_Tab(ab_params, config, folder):
    p   = ab_params['p']
    T0  = ab_params['T0']
    Tab = ab_params['Tab']

    file_params = folder.split("constant_")[-1]

    λd_file = np.loadtxt(f"{folder}/single_optimal_λd-{file_params}.txt")
    λr_file = np.loadtxt(f"{folder}/single_optimal_λr-{file_params}.txt")
    δ_file  = np.loadtxt(f"{folder}/single_optimal_δ-{file_params}.txt")

    # transform p, Tab to indices
    ab_res = int(config['ab_res'])
    T0_max = int(config['T0_max'])
    p_arr  = np.linspace(0, 1,      ab_res)
    T_arr  = np.linspace(0, T0_max, ab_res)

    ip = np.where(abs(p_arr-p)  == np.min(abs(p_arr-p)))
    iT = np.where(abs(T_arr-T0) == np.min(abs(T_arr-T0)))

    # find optimal
    λd_opt = λd_file[ip,iT][0][0]
    λr_opt = λr_file[ip,iT][0][0]
    δ_opt  = δ_file[ip,iT][0][0]

    return λd_opt, λr_opt, δ_opt