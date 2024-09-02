import numpy as np
import matplotlib.pyplot as plt

#folder = args.folder
p   = 0.1
T0  = 0
Tab = 12
T   = T0 + Tab

λd = np.loadtxt(f"data/evolution-T0_{T0:0.0f}-Tab_{Tab:0.0f}/r_arr_specified_competition_average_λd-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")
λr = np.loadtxt(f"data/evolution-T0_{T0:0.0f}-Tab_{Tab:0.0f}/r_arr_specified_competition_average_λr-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")
δ  = np.loadtxt(f"data/evolution-T0_{T0:0.0f}-Tab_{Tab:0.0f}/r_arr_specified_competition_average_δ-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")

# λd = np.loadtxt(f"data/evolution-T0_{T0:0.0f}-Tab_{Tab:0.0f}/competition_average_λd-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")
# λr = np.loadtxt(f"data/evolution-T0_{T0:0.0f}-Tab_{Tab:0.0f}/competition_average_λr-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")
# δ  = np.loadtxt(f"data/evolution-T0_{T0:0.0f}-Tab_{Tab:0.0f}/competition_average_δ-T0_{T0:0.0f}-T_{T:0.0f}-p_{p:0.1f}.txt")


diff_λr = np.diff(λr)
diff_δ  = np.diff(δ)

start = 0
end = 1

fig, ax = plt.subplots(1,3, figsize=(9,2))

ax[0].plot(λd,'-')
ax[1].plot(λr, '-')
ax[2].plot(δ, '-')

for i in [0,1,2]:
    ax[i].set(xscale="log")

# fig, ax = plt.subplots(1,1)
# ax.plot(λr, '-')
# ax.plot(δ, '-')
#ax.plot(diff_λr, '.
# ')
#ax.plot(-diff_λr[diff_λr < 0], '.')
fig.tight_layout()
fig.savefig(f"check_plot_p_{p:0.1f}.png")