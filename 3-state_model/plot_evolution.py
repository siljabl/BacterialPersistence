import sys
import argparse
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append('src')
from config_functions       import read_config
from differential_equations import λ_min, δ_max
from optimal_from_file      import identify_optimal_parameters_const_T0, identify_optimal_parameters_const_Tab

#mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.size"]   = "12"

parser = argparse.ArgumentParser(description='Competition between N species for tot_cycles cycles.')
#parser.add_argument('folder', type=str, help="Folder for heatmap data.")
#parser.add_argument('T0',     type=float, help='application time of antibiotics')
#parser.add_argument('Tab',    type=float, help='duration of antibiotics')
args = parser.parse_args()
tot_cycles = 10_000

p_arr = [0.1, 0.3, 0.5, 0.7, 0.9] #args.p
index = [0, 1, 2, 3, 4]


folder_T0_0 = 'data/constant_T0_0' #args.folder
folder_T0_2 = 'data/constant_Tab_12'
Tab    = 12 #args.Tab

cmap = sns.color_palette('coolwarm', as_cmap=True)
colors = cmap(np.linspace(0,1,5))

sns.set_theme(style='ticks', font_scale=1.1)
fig, ax = plt.subplots(2, 3, figsize=(6.75, 4.5), sharex=True)


for i in index:
    T0 = 2
    T  = T0 + Tab 
    p  = p_arr[i]

    λd = np.loadtxt(f"data/mutation-T0_{T0:0.0f}-Tab_{Tab:0.0f}/competition_average_λd-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")[:tot_cycles]
    λr = np.loadtxt(f"data/mutation-T0_{T0:0.0f}-Tab_{Tab:0.0f}/competition_average_λr-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")[:tot_cycles]
    δ  = np.loadtxt(f"data/mutation-T0_{T0:0.0f}-Tab_{Tab:0.0f}/competition_average_δ-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")[:tot_cycles]

    config = read_config(folder_T0_2)
    ab_params  = {'p': p, 'T0': T0, 'Tab': Tab}
    opt_params = identify_optimal_parameters_const_Tab(ab_params, config, folder_T0_2)

    ax[0,0].plot(λd / T, c=colors[i], alpha=0.7, label=f"{p:0.1}")
    ax[0,1].plot(λr / T, c=colors[i], alpha=0.7, label=f"{p:0.1}")
    ax[0,2].plot(δ,      c=colors[i], alpha=0.7, label=f"{p:0.1}")

    ax[0,0].hlines(opt_params[0] / T, 0, tot_cycles, color=colors[i], ls="dashed")
    ax[0,1].hlines(opt_params[1] / T, 0, tot_cycles, color=colors[i], ls="dashed")
    ax[0,2].hlines(opt_params[2]    , 0, tot_cycles, color=colors[i], ls="dashed")


for i in index:
    T0 = 0
    T  = T0 + Tab 
    p  = p_arr[i]

    λd = np.loadtxt(f"data/mutation-T0_{T0:0.0f}-Tab_{Tab:0.0f}/competition_average_λd-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")[:tot_cycles]
    λr = np.loadtxt(f"data/mutation-T0_{T0:0.0f}-Tab_{Tab:0.0f}/competition_average_λr-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")[:tot_cycles]
    δ  = np.loadtxt(f"data/mutation-T0_{T0:0.0f}-Tab_{Tab:0.0f}/competition_average_δ-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")[:tot_cycles]

    config = read_config(folder_T0_0)
    ab_params  = {'p': p, 'T0': T0, 'Tab': Tab}
    opt_params = identify_optimal_parameters_const_T0(ab_params, config, folder_T0_0)
    
    ax[1,0].plot(λd / T, c=colors[i], alpha=0.7, label=f"{p:0.1}")
    ax[1,1].plot(λr / T, c=colors[i], alpha=0.7, label=f"{p:0.1}")
    ax[1,2].plot(δ,      c=colors[i], alpha=0.7, label=f"{p:0.1}")

    ax[1,0].hlines(opt_params[0] / T, 0, tot_cycles, color=colors[i], ls="dashed")
    ax[1,1].hlines(opt_params[1] / T, 0, tot_cycles, color=colors[i], ls="dashed")
    ax[1,2].hlines(opt_params[2]    , 0, tot_cycles, color=colors[i], ls="dashed")


ax[0,0].set(title=r"$\langle \lambda \rangle ~/~ T$")
ax[0,1].set(title=r"$\langle \omega  \rangle ~/~ T$")
ax[0,2].set(title=r"$\langle \delta  \rangle$")

for j in range(3):
    ax[1,j].set(xlabel=r"Cycle $\times 10^3$")

for i in range(2):
    ax[i,0].set(ylim=[0, 1.1],  yscale="linear", xscale="linear")
    ax[i,1].set(ylim=[0, 1.1],  yscale="linear", xscale="linear")
    ax[i,2].set(ylim=[0, 0.11], xlim=[0,1e4], yscale="linear", xscale="linear")


# Define a custom tick label formatting function
def format_tick_label(x, pos):
    return f"{x/1000:.0f}"

plt.gca().xaxis.set_major_formatter(mpl.ticker.FuncFormatter(format_tick_label))

# fig.tight_layout(rect=[0, 0, 1, 0.85], h_pad=3)
# ax[0,1].legend(loc='upper center',
#                bbox_to_anchor=(.5, 1.9),
#                ncol=5,
#                frameon=False,
#                fancybox=True,
#                handlelength=1,
#                shadow=False,
#                title=r"$p$",
#                alignment='left')

fig.tight_layout(rect=[0, 0, 0.9, 1], h_pad=3)
ax[0,1].legend(loc='upper center',
               bbox_to_anchor=(3.2, 0.5),
               ncol=1,
               reverse=True,
               frameon=False,
               fancybox=True,
               handlelength=1,
               title=r"$p$",
               alignment='center')

sns.despine()
fig.savefig(f"figs/mutation_average/mutation_average_parameters.png")