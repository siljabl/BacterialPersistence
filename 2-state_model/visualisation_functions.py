# Two functions for plotting data

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from model_equations import ode_grow, ode_kill
from model_equations import analytical_population, analytical_growth, analytical_death
from model_equations import n0, S0, f, n_min
from simulation_functions import fitness, winner_params, optimal_competitor_analytical, optimal_competitor_stochastic


font = {'family': 'Times New Roman',
        'weight': 'normal',
        'size': 20}

mpl.rc('font', **font)
fs_text = 14


###################
## Visualisation ##
###################
# plotting behavior of system, both with ode and analytical
def plot_single_sample(bac_args, ab_args, sim_args, r_arr, winner=True, compare=False, save_fig=False):
    # importing
    if winner:
        bac_args = winner_params(bac_args, ab_args)
    p, T0, Tab = ab_args[0:3]
    tot_cycles = sim_args[3]

    NA = np.newaxis
    T, T_cycle = T0 + Tab, 40

    # Empty arrays to for plot
    n0 = 1
    n0_d, nt_d = np.array([n0, n0, 0, 0, S0]), [[], [], [], [], []]  # arrays for ode solver
    n0_a, nt_a = [np.array([n0, n0]), np.zeros(2), S0], [[], [], [], []]  # arrays for analytic solution
    ode_time, a_time = [], []
    full_a, full_t = [[], []], []

    # solving dynamics
    for ic in range(tot_cycles):
        r = r_arr[ic]

        # with antibiotics
        if r < p:
            # solving ODEs
            sol_cycle1 = solve_ivp(ode_grow, [0, T0], n0_d, args=bac_args[0:2], max_step=0.1)
            sol_cycle2 = solve_ivp(ode_kill, [T0, T], sol_cycle1.y[:, -1], args=bac_args[0:2], max_step=0.1)
            if threshold:
                dead = (sol_cycle2.y[0:4] < n_min)
                sol_cycle2.y[0:4][dead] = 0
            sol_cycle3 = solve_ivp(ode_grow, [T, T_cycle], sol_cycle2.y[:, -1], args=bac_args[0:2], max_step=0.1)

            # assembling arrays
            t_temp = np.concatenate([sol_cycle1.t, sol_cycle2.t, sol_cycle3.t])
            n_temp = np.concatenate([sol_cycle1.y.T, sol_cycle2.y.T, sol_cycle3.y.T]).T

            # analytical solution
            # analytic1 = analytical_growth(sol_cycle1.t[:, NA], n0_a, *bac_args[2:4])
            # a_T0 = [analytic1[i][-1] for i in range(3)]
            # analytic2 = analytical_death(sol_cycle2.t[:, NA] - T0, a_T0, *bac_args[2:6])
            # a_T = [analytic2[i][-1] for i in range(3)]
            if threshold:
                for i in range(2):
                    dead = (a_T[i] < n_min)
                    a_T[i][dead] = 0
            # analytic3 = analytical_growth(sol_cycle3.t[:, NA][sol_cycle3.y[4] > 0] - T, a_T, *bac_args[2:4])

            t = max(sol_cycle3.t[sol_cycle3.y[4] > 0])
            full_t.append(t + ic * T_cycle)
            # full_a[0].append(analytical_population(t, T, T0, n0_a[0][0], bac_args[0][0], bac_args[1][0]))
            # full_a[1].append(analytical_population(t, T, T0, n0_a[0][1], bac_args[0][1], bac_args[1][1]))
            # assembling arrays
            # for i in range(2):
                # nt_a[2 * i].append(np.concatenate([analytic1[i][:, 0], analytic2[i][:, 0], analytic3[i][:, 0], [analytic3[i][-1, 0]]]))
                # nt_a[2 * i + 1].append(np.concatenate([analytic1[i][:, 1], analytic2[i][:, 1], analytic3[i][:, 1], [analytic3[i][-1, 1]]]))

        # without antibiotics
        else:
            # solving ODE
            sol_cycle = solve_ivp(ode_grow, [0, T_cycle], n0_d, args=bac_args[0:2], max_step=0.1)
            t_temp, n_temp = sol_cycle.t, sol_cycle.y

            t = max(sol_cycle.t[sol_cycle.y[4] > 0])
            full_t.append(t + ic * T_cycle)
            # full_a[0].append(analytical_population(t, 0, 0, n0_a[0][0], bac_args[0][0], bac_args[1][0]))
            # full_a[1].append(analytical_population(t, 0, 0, n0_a[0][1], bac_args[0][1], bac_args[1][1]))

            # analytical solution
            # analytical = analytical_growth(t_temp[:, NA][n_temp[4] > 0], n0_a, *bac_args[2:4])
            # for i in range(2):
                # nt_a[2 * i].append(np.concatenate([analytical[i][:, 0], [analytical[i][-1, 0]]]))
                # nt_a[2 * i + 1].append(np.concatenate([analytical[i][:, 1], [analytical[i][-1, 1]]]))

        # saving data
        ode_time.append(t_temp + ic * T_cycle)
        a_time.append(ode_time[ic][n_temp[4] > 0])
        a_time.append([(ic + 1) * T_cycle])
        for i in range(5):
            nt_d[i].append(n_temp[i])

        # preparing new cycle
        n0_d = f * np.array([nt_d[0][ic][-1] + nt_d[2][ic][-1], nt_d[1][ic][-1] + nt_d[3][ic][-1], 0, 0, S0 / f])
        # n0_a = [f * np.array([nt_a[0][ic][-1] + nt_a[2][ic][-1], nt_a[1][ic][-1] + nt_a[3][ic][-1]]), np.zeros(2), S0]
        if threshold:
            dead = (n0_d[0:4] < n_min)
            n0_d[0:4][dead] = 0
            # for i in range(2):
                # dead = (n0_a[i] < n_min)
                # n0_a[i][dead] = 0

    # concatenating arrays for plotting
    for i in range(4):
        nt_d[i] = np.concatenate(nt_d[i])
        # nt_a[i] = np.concatenate(nt_a[i])
    S = np.concatenate(nt_d[4])
    ode_time = np.concatenate(ode_time)
    # a_time = np.concatenate(a_time)

    if winner: labels = ['winner', 'competitor']
    else: labels = ['(3.6, 0)', '(0.01, 0)']
    species = ['(' + str(bac_args[0][i].round(1)) + ', ' + str(bac_args[1][i].round(2)) + ')' for i in range(2)]
    color = ['tab:orange', 'tab:green', 'tab:blue']
    title = "p = " + str(np.round(p, 2)) + ",  Tab = " + str(np.round(Tab, 2)) + ",  " + labels[0] + ": " + species[0] + ",  " + labels[1] + ": " + species[1]

    if compare:
        fig, ax = plt.subplots(2, 1, figsize=(8, 3.5), sharex="all", sharey="all")
        fig.suptitle(title)

        for i in range(2):
            ax[i].plot(ode_time, S, alpha=0.4, label="substrate")
        for i in range(2):
            ax[i].plot(ode_time, nt_d[i] + nt_d[i+2], color=color[0], label='ODE')
            ax[i].plot(a_time, nt_a[i] + nt_a[i+2], '--', color=color[1], label='analytic parts')
            ax[i].plot(full_t, full_a[i], 'x', color=color[2], label="full analytic")

        ax[0].set(ylabel="Log(population)", title=labels[0])
        ax[1].set(xlabel="Time [h]", ylabel="Log(population)", title=labels[1])
        ax[0].set(yscale="log", ylim=(1, 10 ** 10))
        ax[0].legend(bbox_to_anchor=(1, 1), loc="best")
        ax[1].legend(bbox_to_anchor=(1, 1), loc="best")

    else:
        fig, ax = plt.subplots(1, 1, figsize=(12, 4), sharex="all", sharey="all")
        ax.plot(ode_time, S, lw=3, alpha=0.5, label="substrate")
        for i in range(2):
            ax.plot(ode_time, nt_d[i] + nt_d[i + 2], color=color[i], lw=3, label=r"($\lambda$, $\delta$) = " + labels[i])
            # ax.plot(a_time, nt_a[i+2], color=color[i], label=labels[i])
        ax.set(xlabel="Time [h]", ylabel="Log(population)")#, title=title)
        ax.set(yscale="log", ylim=(10**(-2), 10 ** 10))
        fig.tight_layout(rect=[0, 0, 1, 0.9])
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.4),
              ncol=3, fancybox=True, shadow=False)


    #fig.tight_layout()
    fig.show()
    if save_fig:
        fig.savefig("../../../figures/model2/ode_example.png", dpi=100)



# plotting extinction frequency from competition data
def extinction_frequency(ab_args, indices, save_fig=False):
    cmap = mpl.cm.get_cmap('cividis')
    # importing data
    ip, it = indices
    p, T0, Tab = ab_args
    p = np.round(p, 1)
    Tab = np.round(Tab, 0)
    ab_label = [str(np.round(ab_arg, 1)) for ab_arg in ab_args[0:3]]
    lag_opt = np.loadtxt('../../../data/model2/low_resolution/optimal_lag-T0' + str(T0))
    del_opt = np.loadtxt('../../../data/model2/low_resolution/optimal_delta-T0' + str(T0))
    w_params = '(' + str(np.round(lag_opt[ip+1, it], 2)) + ', ' + str(np.round(del_opt[ip+1, it], 2)) + ')'
    
    parameters = "p" + str(int(10*p)*10) + "-T0" + str(int(T0)) + "-Tab" + str(int(Tab))
    p_ext_winner = np.loadtxt("../../../data/model2/competition_two_species/extinction_frequency/optimal_extinction_prob-" + parameters)
    p_ext_competitor = np.loadtxt("../../../data/model2/competition_two_species/extinction_frequency/competitor_extinction_prob-" + parameters)
    lag_opt = np.loadtxt('../../../data/model2/competition_two_species/competition_lag-T0'+str(T0))
    del_opt = np.loadtxt('../../../data/model2/competition_two_species/competition_delta-T0'+str(T0))
    optimal_comp = '(' + str(np.round(lag_opt[ip, it], 2)) + ', ' + str(np.round(del_opt[ip, it], 2)) + ')'
    # optimal_comp = '(?, ?)'

    # setting up figure
    fig, ax = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    fig.suptitle(r'Extinction frequency, $p$ = ' + ab_label[0] + r', $T_0$ = ' + ab_label[1] + r', $T_{ab}$ = ' + ab_label[2])
    fig.suptitle('$p$ = ' + ab_label[0] + r', $T_0$ = ' + ab_label[1] + r', $T_{ab}$ = ' + ab_label[2])
    ax[0].set(xlabel=r'$\lambda_{II}$ [h]', ylabel=r'$\delta_{II}$', title=r'$(\lambda, \delta)^*_I$ = ' + w_params)
    ax[1].set(xlabel=r'$\lambda_{II}$ [h]', title=r'$(\lambda, \delta)^*_{II}$ = ' + optimal_comp)

    # plotting
    im = ax[0].imshow(p_ext_winner, vmin=0, vmax=1, aspect='auto', cmap=cmap, origin='lower', extent=[0.01, T0 + 12, 0, 0.1])
    im = ax[1].imshow(p_ext_competitor, vmin=0, vmax=1, aspect='auto', cmap=cmap, origin='lower', extent=[0.01, T0 + 12, 0, 0.1])
    # ax[0].plot(lag_opt[ip, it], del_opt[ip, it], 's', lw=3, color='r', markersize=8)
    # ax[1].plot(lag_opt[ip, it], del_opt[ip, it], 's', lw=3, color='r', markersize=8)

    # saving
    fig.tight_layout(rect=[0, 0, 1, 1.02])
    fig.colorbar(im, ax=ax[0:2])
    fig.show()
    if save_fig:
        fig.savefig("../../../figures/model2/extinction_frequency_p" + ab_label[0] + "_T0" + ab_label[1] + "_T" + ab_label[2] + ".png", dpi=100)

    fig, ax = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    # ax[0].plot(np.linspace(0, T0 + 12, 100), p_ext_winner[0], '.-', label=r'winner, $\delta$ = ' + str(0.1*0))
    # ax[0].plot(np.linspace(0, T0 + 12, 100), p_ext_competitor[0], '.-', label=r'competitor, $\delta$ = ' + str(0.1 * 0))
    # ax[0].plot(np.linspace(0, T0 + 12, 100), p_ext_winner[60], '.-', label=r'winner, $\delta$ = ' + str(0.1 * 0.6))
    # ax[0].plot(np.linspace(0, T0 + 12, 100), p_ext_competitor[60], '.-', label=r'competitor, $\delta$ = ' + str(0.1 * 0.6))
    #
    # ax[1].plot(np.linspace(0, 0.25, 100), p_ext_winner[:, 20], '.-', label=r'winner, $\lambda$ = ' + str(np.round((T0+Tab) * 0.2, 1)))
    # ax[1].plot(np.linspace(0, 0.25, 100), p_ext_competitor[:, 20], '.-', label=r'competitor, $\lambda$ = ' + str(np.round((T0+Tab) * 0.2, 1)))
    # ax[1].plot(np.linspace(0, 0.25, 100), p_ext_winner[:, 60], '.-', label=r'winner, $\lambda$ = ' + str((T0+Tab) * 0.6))
    # ax[1].plot(np.linspace(0, 0.25, 100), p_ext_competitor[:, 60], '.-', label=r'competitor, $\lambda$ = ' + str((T0+Tab) * 0.6))
    ax[0].plot(np.linspace(0.01, T0 + 12, 100), p_ext_competitor[20], '.-', label=r'comp, $\delta$ = ' + str(0.01 * 20))
    ax[0].plot(np.linspace(0.01, T0 + 12, 100), p_ext_winner[20], '.-', label=r'winner, $\delta$ = ' + str(0.01 * 20))

    ax[0].set(xlabel=r'$\lambda$ [h]', ylabel='Extinction frequency')
    ax[1].set(xlabel=r'$\delta$')
    ax[0].legend()
    ax[1].legend()
    fig.tight_layout()
    fig.show()


