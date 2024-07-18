import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_bacterial_parameters(bac_params, folder):
    fig, ax = plt.subplots(1,3, figsize=(9,3))
    ax[0].plot(bac_params['λd'], '.-')
    ax[1].plot(bac_params['λr'], '.-')
    ax[2].plot(bac_params['δ'],  '.-')

    ax[0].set(title="λd")
    ax[1].set(title="λr")
    ax[2].set(title="δ")

    fig.tight_layout()
    fig.savefig(f"{folder}/bacterial_parameters.png")


def plot_cycles(sol_cycles, bac_params, ab_params, sim_params, folder):
    species = sol_cycles[0]
    substrate = sol_cycles[1]
    time = sol_cycles[2]

    λd = bac_params['λd']
    λr = bac_params['λr']
    δ  = bac_params['δ']

    max_cycle = sim_params['tot_cycles']

    cycles = np.arange(int(max_cycle / 10)-1, max_cycle, int(max_cycle / 10))
    colors = mpl.cm.jet(np.linspace(0,1,20))

    fig, axes = plt.subplots(2, 5, figsize=(10,4), sharex=True, sharey=True)
    ax = axes.flatten()

    ic = 0
    for c, ax_c in zip(cycles, ax):
        # plotting amount of nutrients as blue background
        ax_c.fill_between(time[ic], 0, substrate[ic], color="lightskyblue")

        # looping through all species
        for i in range(len(species[ic])-1):
            
            # only plotting species that are alive at start of cycle
            if species[ic][i][0] > 1:
                ax_c.plot(time[ic], species[ic][i], color=colors[i%12], label=f"({λd[i]:0.2f}, {λr[i]:0.2f}, {δ[i]:0.2f})")

        ax_c.set(title=f"cycle #{c}", yscale="linear")
        ic += 1
    
    fig.tight_layout(rect=(0,0,1,0.95))    
    # ax[2].legend(loc='upper center', bbox_to_anchor=(0.5, 1.5),
    #               ncol=5, fancybox=True, shadow=True)

    T0 = ab_params['T0']
    T  = ab_params['T']
    p  = ab_params['p']

    fig.savefig(f"{folder}/feast-famine_cycles-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.png")    