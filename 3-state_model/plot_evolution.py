import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('src')
from differential_equations import λ_min, δ_max

parser = argparse.ArgumentParser(description='Competition between N species for tot_cycles cycles.')
# parser.add_argument('p',   type=float, help='frequency of antibiotics.')
parser.add_argument('T0',  type=float, help='application time of antibiotics')
parser.add_argument('Tab', type=float, help='duration of antibiotics')
args = parser.parse_args()

p_arr = [0.1] #args.p
T0  = args.T0
Tab = args.Tab
T   = T0 + Tab

# fig_pop,   ax_pop   = plt.subplots(1, 3, figsize=(11,3))
fig_param, ax_param = plt.subplots(1, 3, figsize=(11,3))


for p in p_arr:
    λd = np.loadtxt(f"data/competition_average/average_λd-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")
    λr = np.loadtxt(f"data/competition_average/average_λr-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")
    δ  = np.loadtxt(f"data/competition_average/average_δ-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")

    # p_dominant = np.loadtxt(f"data/competition_average/dominant_species-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")
    # n_extinct  = np.loadtxt(f"data/competition_average/number_of_extinctions-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")
    # p_dists = np.loadtxt(f"data/competition_average/population_distributions-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")

    ########################
    ## Average parameters ##
    ########################
    # fig_param, ax_param = plt.subplots(1, 3, figsize=(11,3))

    ax_param[0].plot(λd / T, label=f"p={p:0.1}")
    ax_param[1].plot(λr / T, label=f"p={p:0.1}")
    ax_param[2].plot(δ, label=f"p={p:0.1}")


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

ax_param[0].set(xlabel="Cycle number", ylabel="λd / T")
ax_param[1].set(xlabel="Cycle number", ylabel="λr / T")
ax_param[2].set(xlabel="Cycle number", ylabel="δ")

ax_param[0].set(ylim=[0, 1.1],     yscale="linear")
ax_param[1].set(ylim=[0, 1.1],     yscale="linear")
ax_param[2].set(ylim=[0, δ_max], yscale="linear")

fig_param.tight_layout()
fig_param.savefig(f"figs/competition_average/average_parameters-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.png")

# fig_pop.tight_layout()
# fig_pop.savefig(f"figs/competition_average/populations-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.png")
