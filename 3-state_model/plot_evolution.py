import sys
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

sys.path.append('src')
from config_functions       import read_config
from differential_equations import λ_min, δ_max
from optimal_from_file      import identify_optimal_parameters_const_T0, identify_optimal_parameters_const_Tab

parser = argparse.ArgumentParser(description='Competition between N species for tot_cycles cycles.')
parser.add_argument('folder',      type=str, help="Folder for heatmap data.")
parser.add_argument('T0',  type=float, help='application time of antibiotics')
parser.add_argument('Tab', type=float, help='duration of antibiotics')
args = parser.parse_args()

p_arr  = [0.1, 0.3, 0.5, 0.7, 0.9] #args.p
#p_arr  = [0.9, 0.7, 0.5, 0.3, 0.1] #args.p

folder = args.folder
T0     = args.T0
Tab    = args.Tab
T      = T0 + Tab

tot_cycles = 10_000
config = read_config(folder)


colors = mpl.cm.jet(np.linspace(0,1,5))
# fig_pop,   ax_pop   = plt.subplots(1, 3, figsize=(11,3))
fig_param, ax_param = plt.subplots(1, 3, figsize=(6.75, 2.75))

index = [0, 1, 2, 3, 4]
for i in index:
    p = p_arr[i]
    λd = np.loadtxt(f"data/mutation-T0_{T0:0.0f}-Tab_{Tab:0.0f}/competition_average_λd-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")[:tot_cycles]
    λr = np.loadtxt(f"data/mutation-T0_{T0:0.0f}-Tab_{Tab:0.0f}/competition_average_λr-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")[:tot_cycles]
    δ  = np.loadtxt(f"data/mutation-T0_{T0:0.0f}-Tab_{Tab:0.0f}/competition_average_δ-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")[:tot_cycles]

    ab_params  = {'p': p, 'T0': T0, 'Tab': Tab}
    opt_params = identify_optimal_parameters_const_T0(ab_params, config, folder)
    #opt_params = identify_optimal_parameters_const_Tab(ab_params, config, folder)

    # p_dominant = np.loadtxt(f"data/competition_average/dominant_species-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")
    # n_extinct  = np.loadtxt(f"data/competition_average/number_of_extinctions-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")
    # p_dists = np.loadtxt(f"data/competition_average/population_distributions-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")

    ########################
    ## Average parameters ##
    ########################
    # fig_param, ax_param = plt.subplots(1, 3, figsize=(11,3))

    ax_param[0].plot(λd / T, c=colors[i], alpha=0.9, label=f"p={p:0.1}")
    ax_param[1].plot(λr / T, c=colors[i], alpha=0.9, label=f"p={p:0.1}")
    ax_param[2].plot(δ,      c=colors[i], alpha=0.9, label=f"p={p:0.1}")

    ax_param[0].hlines(opt_params[0] / T, 0, tot_cycles, color=colors[i], ls="dashed")
    ax_param[1].hlines(opt_params[1] / T, 0, tot_cycles, color=colors[i], ls="dashed")
    ax_param[2].hlines(opt_params[2]    , 0, tot_cycles, color=colors[i], ls="dashed")


    ###########################
    ## Antibiotic parameters ##
    ###########################
    # fig_pop, ax_pop = plt.subplots(1, 3, figsize=(11,3))
    # ax_pop[0].plot(p_dominant, label=f"p={p:0.1}")
    # for i in range(len(p_dists[0])):
    #     ax_pop[1].plot(p_dists[:,i], label=f"{i}")
    # ax_pop[2].plot(n_extinct,  label=f"p={p:0.1}")

    #ax_pop[1].set(xscale="log", yscale="log", ylim=(0.01, 2))

# ax_pop[1].legend()
# for i in range(3):
#     ax_param[i].legend()
#     ax_pop[i].legend()

ax_param[0].set(xlabel="Cycle number", title=r"$\langle \lambda \rangle ~/~ T$")
ax_param[1].set(xlabel="Cycle number", title=r"$\langle \omega \rangle ~/~ T$")
ax_param[2].set(xlabel="Cycle number", title=r"$\langle \delta \rangle$")

ax_param[0].set(ylim=[0, 1.1], yscale="linear", xscale="linear")
ax_param[1].set(ylim=[0, 1.1], yscale="linear", xscale="linear")
ax_param[2].set(ylim=[0, 0.1], yscale="linear", xscale="linear")

fig_param.tight_layout(rect=[0, 0, 1, 0.85])
ax_param[1].legend(loc='upper center', bbox_to_anchor=(.5, 1.5),
          ncol=5, fancybox=True, shadow=False)


fig_param.savefig(f"figs/mutation_average/mutation_average_parameters-T0_{T0:0.0f}-T_{T:0.0f}.png")

# fig_pop.tight_layout()
# fig_pop.savefig(f"figs/competition_average/populations-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.png")
