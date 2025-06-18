import sys
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append('src')
#from config_functions  import read_config
from optimal_from_file import optimal_parameters_from_data

tot_cycles = 10_000
p_arr = [0.1, 0.3, 0.5, 0.7, 0.9]
index = [2]

dir = "data/mutation/mutation_rate-1e-03"
folder_T0_0 = 'data/constant_T0_0'
folder_T0_2 = 'data/constant_Tab_12'


# Configure figure
cmap = sns.color_palette('coolwarm', as_cmap=True)
colors = cmap(np.linspace(0,1,5))

sns.set_theme(style='ticks', font_scale=1.2)
fig, ax = plt.subplots(2, 2, figsize=(6.75, 5), sharex=True)

for i in index:
    p  = p_arr[i]
    Tab = 12
    T0  = 5
    T  = T0 + Tab

    λ = np.loadtxt(f"{dir}/average_lag-p{p*100:0.0f}-T0{T0:0.0f}-T{T:0.0f}-min")
    δ = np.loadtxt(f"{dir}/average_delta-p{p*100:0.0f}-T0{T0:0.0f}-T{T:0.0f}-min")  

    ax[0,0].plot(λ / T, c=colors[i], alpha=0.7, label=f"{p:0.1}")
    ax[1,0].plot(δ,     c=colors[i], alpha=0.7, label=f"{p:0.1}")

    Tab = 18
    T0  = 5
    T  = T0 + Tab

    λ = np.loadtxt(f"{dir}/average_lag-p{p*100:0.0f}-T0{T0:0.0f}-T{T:0.0f}-min")
    δ = np.loadtxt(f"{dir}/average_delta-p{p*100:0.0f}-T0{T0:0.0f}-T{T:0.0f}-min")

    ax[0,1].plot(λ / T, c=colors[i], alpha=0.7, label=f"{p:0.1}")
    ax[1,1].plot(δ,     c=colors[i], alpha=0.7, label=f"{p:0.1}")




ax[0,0].set(ylabel=r"$\langle \lambda \rangle ~/~ T$")
ax[1,0].set(ylabel=r"$\langle \delta  \rangle$")

for i in range(2):
    ax[1,i].set(xlabel=r"Cycle $\times 10^3$")
    ax[0,i].set(ylim=[0, 1.1],  yscale="linear", xscale="linear")
    ax[1,i].set(ylim=[0, 0.2], xlim=[0,1e4], yscale="linear", xscale="linear")


def format_tick_label(x, pos):
    return f"{x/1000:.0f}"

plt.gca().xaxis.set_major_formatter(mpl.ticker.FuncFormatter(format_tick_label))


fig.tight_layout(rect=[0, 0, 0.87, 1], h_pad=3)
ax[1,1].legend(loc='upper center',
               bbox_to_anchor=(1.3, 1.75),
               ncol=1,
               reverse=True,
               frameon=False,
               fancybox=True,
               handlelength=1,
               title=r"$p$",
               alignment='center')

sns.despine()
fig.savefig("figs/mutation/T0_5_Tab_12-18.png", dpi=300)