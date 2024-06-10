import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

font = {'family': 'Times New Roman',
        'weight': 'normal',
        'size': 20}
mpl.rc('font', **font)

folder = "/home/silja/Documents/thesis/thesis_git/type_II/data/competition_plot/comp_"
# folder = "/home/silja/Documents/thesis/thesis_git/type_II/data/single_opt/single_"


T_labels = ['T0', 'Tab']
path = folder
ab_res = 100

T0 = 0.0
Tab = np.array([6.0])
Nt = len(Tab)
t_index = np.array([int(ab_res * t / 24) for t in Tab])

cases = ['I', 'II', 'III']
mut_seeds = ['first', 'optimal', 'max']
models = [1, 2, 3]


#######################
## Folder parameters ##
#######################
mut_seed = mut_seeds[0]
case = cases[1]
model = models[1]

data_path = "model_" + str(model) + "/case_" + case + "/data/server/seed_" + mut_seed
fig_path = "model_" + str(model) + "/case_" + case + "/figs/seed_" + mut_seed


####################
## Importing data ##
####################
T = (T0 + Tab[0])

lag_avrg, delta_avrg = [], []
p_range = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
for pi in p_range:
    lag_avrg.append(np.loadtxt(data_path + "-average_lag-p" + str(pi) + "-T0" + str(int(T0)) + "-T" + str(int(T))))
    delta_avrg.append(np.loadtxt(data_path + "-average_delta-p" + str(pi) + "-T0" + str(int(T0)) + "-T" + str(int(T))))

lag_p = np.zeros_like(p_range)
lag_std = np.zeros_like(p_range)
delta_p = np.zeros_like(p_range)
delta_std = np.zeros([2, len(p_range)])

for i in range(len(p_range)):
    lag_p[i] = np.median(lag_avrg[i][1000:])
    lag_std[i] = np.std(lag_avrg[i][1000:])
    delta_p[i] = np.median(delta_avrg[i][1000:])
    delta_std[0, i] = np.std(delta_avrg[i][1000:])
    delta_std[1, i] = np.std(delta_avrg[i][1000:])
    if delta_p[i] - delta_std[0,i] < 0:
        delta_std[0,i] = delta_p[i]


#######################
## PREPARING TO PLOT ##
#######################
lag_opt_p = np.zeros([Nt, ab_res])
del_opt_p = np.zeros([Nt, ab_res])

i = 0
for t in Tab:
    lag_opt_p[i] = np.loadtxt(path + 'lag-T0' + str((T0)))[:, t_index[i % Nt]]
    del_opt_p[i] = np.loadtxt(path + 'delta-T0' + str((T0)))[:, t_index[i % Nt]]
    i += 1

p_arr = np.linspace(0, 1, ab_res)
lc = ['firebrick', 'steelblue']
le = ['firebrick', 'steelblue']
c = 1
j = 0

fig, ax = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
ax[0].set(xlabel=r"$p$", ylabel=r"$\langle \lambda \rangle_c$")
ax[1].set(xlabel=r"$p$", ylabel=r"$\langle \delta \rangle_c$")

for i in range(Nt):
    mask0 = (lag_opt_p[i] < 0.1)
    mask1 = (lag_opt_p[i] > 0.1)

    ax[0].plot(p_arr[mask0], lag_opt_p[i][mask0], '-', fillstyle='none', color=lc[c], alpha=0.4, lw=3, markeredgewidth=3)
    ax[0].plot(p_arr[mask1][0:-1:4], lag_opt_p[i][mask1][0:-1:4], '-', fillstyle='none', color=lc[c], alpha=0.4, lw=3, markeredgewidth=3)
    ax[0].errorbar(p_range, lag_p, yerr=lag_std, fmt='o', fillstyle='full', color=le[c], elinewidth=1, capsize=3, label=r"$T_0$ = " + str(T0) + r",  $T_{AB}$ = " + str(T-T0))

    ax[1].plot(p_arr[mask0], del_opt_p[i][mask0], '-', color=lc[c], alpha=0.5, lw=3)
    ax[1].plot(p_arr[mask1][0:-1:4], del_opt_p[i][mask1][0:-1:4], '-', color=lc[c], alpha=0.5, lw=3)
    ax[1].errorbar(p_range, delta_p, yerr=delta_std, fmt='o', fillstyle='full', color=le[c], elinewidth=1, capsize=3, label=r"$T_0$ = " + str(T0) + r",  $T_{AB}$ = " + str(T-T0))

fig.tight_layout(rect=[0, 0, 1, 0.85])
ax[1].legend(loc='upper center', bbox_to_anchor=(-0.18, 1.3), ncol=2, fancybox=True, shadow=False)
fig.show()
#
fig.savefig("model_2/case_II/figs/average_p_ext.pdf")

#fig1 = fig

