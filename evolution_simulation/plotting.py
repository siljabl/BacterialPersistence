import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

folder = ""
folder_fig = "/"

cases = ['I', 'II', 'III']
mut_seeds = ['first', 'optimal', 'max']
models = [1, 2, 3]

# setting parameters..
font = {'family': 'Times New Roman',
        'weight': 'normal',
        'size': 20}
# plt.rcParams.update({"text.usetex": True})
mpl.rc('font', **font)


#######################
## Folder parameters ##
#######################
ifig = 2
case = cases[2]
mut_seed = mut_seeds[ifig]
model = models[1]
save_fig = True

data_path = "model_" + str(model) + "/case_" + case + "/data" + folder + "/seed_" + mut_seed
fig_path = "model_" + str(model) + "/case_" + case + "/figs"+folder_fig+"seed_" + mut_seed

cmap = mpl.cm.get_cmap('CMRmap')


####################
## Importing data ##
####################
T0 = 5.0
T = T0 + 12
p_range = [0.1, 0.3, 0.5, 0.7, 0.9]
colors = [cmap(i/len(p_range)) for i in range(len(p_range))]

lag_avrg, lag_dom = [], []
delta_avrg, delta_dom = [], []
omega_avrg, omega_dom = [], []
extinctions = []
for pi in p_range:
    lag_avrg.append(np.loadtxt(data_path + "-average_lag-p" + str(pi) + "-T0" + str(T0) + "-T" + str(T)))
    lag_dom.append(np.loadtxt(data_path + "-dominant_lag-p" + str(pi) + "-T0" + str(T0) + "-T" + str(T)))

    omega_avrg.append(np.loadtxt(data_path + "-average_omega-p" + str(pi) + "-T0" + str(T0) + "-T" + str(T)))
    omega_dom.append(np.loadtxt(data_path + "-dominant_omega-p" + str(pi) + "-T0" + str(T0) + "-T" + str(T)))


    if case == "III":
        extinctions.append(np.loadtxt(data_path + "-extinctions-p" + str(pi) + "-T0" + str(T0) + "-T" + str(T)))

# p_crit = [0.2, 0.3, 0.4]
# c_crit = ['firebrick', 'orangered', 'gold']
# lag_crit = []
# for pi in p_crit:
#     lag_crit.append(np.loadtxt(data_path + "-average_lag-p" + str(pi) + "-T0" + str(T0) + "-T" + str(T)))


##############
## Plotting ##
##############
title = "Model " + str(model) + " - Case " + case
colors = ['purple', 'red', 'orange', 'yellowgreen', 'green', 'cyan', 'blue', 'black', 'magenta']
# colors = ['peachpuff', 'indigo', 'crimson', 'seagreen', 'blue', 'black']

tot_cycles = len(lag_avrg[0])
cycles = np.arange(1, tot_cycles+1)
xTicks = ['0', r'$10^3$']
# plotting average lag and delta
fig, ax = plt.subplots(1, 3, figsize=(12, 4.5))
ax[0].set(xlabel='Cycle number', ylabel=r'$\langle\lambda\rangle$') #, title=title)
ax[2].set(xlabel='Cycle number', ylabel=r'$\langle\delta\rangle$') #, title=title)
ax[1].set(xlabel='Cycle number', ylabel=r'$\langle\omega\rangle$') #, title=title)

i = 0
for pi in p_range:
    ax[0].plot(cycles, lag_avrg[i], color=colors[i], label="p = "+str(pi))
    ax[2].plot(cycles, delta_avrg[i], color=colors[i], label="p = "+str(pi))
    ax[1].plot(cycles, omega_avrg[i], color=colors[i], label="p = "+str(pi))
    i += 1

ax[0].set(yscale='log', ylim=(0.01, 20))
ax[1].set(yscale='log', ylim=(0.01, 20))
ax[2].set(yscale='log', ylim=(5*10**(-4), 0.2))
# ax[0].legend()
# fig.tight_layout()
fig.tight_layout(rect=[0, 0, 1, 0.8])
fig.subplots_adjust(wspace=0.4, left=0.08, right=0.98)
ax[1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.5),
                ncol=5, fancybox=True, shadow=False)
fig.show()
if save_fig:
    fig.savefig(fig_path + "-averages-T0" + str(int(T0)) + "-T" + str(int(T)))

# plotting dominant lag and delta
fig, ax = plt.subplots(1, 3, figsize=(12, 4))
ax[0].set(xlabel='Cycle number', ylabel='Dominant lag time [h]', title=title)
ax[1].set(xlabel='Cycle number', ylabel='Dominant delta', title=title)
ax[2].set(xlabel='Cycle number', ylabel='Dominant omega [h]', title=title)

i = 0
for pi in p_range:
    ax[0].plot(cycles, lag_dom[i],   color=colors[i], label="p = "+str(pi))
    ax[1].plot(cycles, delta_dom[i], color=colors[i], label="p = "+str(pi))
    ax[2].plot(cycles, omega_dom[i], color=colors[i], label="p = "+str(pi))
    i += 1

ax[0].set(yscale='log', ylim=(0.009, 10))
ax[1].set(yscale='linear')
ax[0].legend()
fig.tight_layout()
fig.show()
if save_fig:
    fig.savefig(fig_path + "-dominants-T0" + str(T0) + "-T" + str(T))



if case == 'III':
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    ax.set(xlabel='Cycle number', ylabel='# Extinctions', title=title)

    i = 0
    for pi in p_range:
        ax.plot(cycles, extinctions[i], color=colors[i], label="p = " + str(pi))
        i += 1

    ax.set(xscale='log', yscale='linear')
    ax.legend()
    fig.tight_layout()
    #fig.show()
    if save_fig:
        fig.savefig(fig_path + "-extinctions-T0" + str(T0) + "-T" + str(T))

