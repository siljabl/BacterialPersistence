import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


p_arr = np.linspace(0.1, 0.9, 9)
data_T0_0 = np.loadtxt(f"../data/supplementary/strategy_fitness_T00.0_Tab12.0.txt")
data_T0_2 = np.loadtxt(f"../data/supplementary/strategy_fitness_T02.0_Tab12.0.txt")


sns.set_theme(style='ticks', font_scale=1.2)
fig, ax = plt.subplots(1, 2, figsize=(8, 3.5), sharey=True, sharex=True)

ax[0].plot(p_arr, data_T0_0[:,0], 'v-', label="Triggered")
ax[0].plot(p_arr, data_T0_0[:,1], '^-', label="Spontaneous")
ax[0].set(xlabel=r'$p$',  ylabel=r"$F(\lambda, \omega, \delta; p, T_0, T_{ab})$")

ax[1].plot(p_arr, data_T0_2[:,0], 'v-', label="Triggered")
ax[1].plot(p_arr, data_T0_2[:,1], '^-', label="Spontaneous")
ax[1].set(xlabel=r'$p$')

fig.tight_layout()
ax[0].legend(loc='best', frameon=False)
ax[1].legend(loc='best', frameon=False)

fig.savefig(f"../figs/supplementary/strategy_fitness.png")