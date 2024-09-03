# Plotting lineplots from data
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

font = {'family': 'Times New Roman',
        'weight': 'normal',
        'size': 20}
mpl.rc('font', **font)

T_labels = ['T0', 'Tab']
data_path = '../../../data/model2/'
ab_res = 400
saving = False

##############
## SETTINGS ##
##############
time_label = 'T0'						# time variable that is kept constant. 'Tab' or 'T0'
time_value = 5 							# value of constant variable
time = time_label + str(time_value)

x_var = 'p'							# variable that parameters are plotted against. 'p', 'Tab' or 'T0'
line_var = 'Tab'						# variable that is kept constant for each line. 'p', 'Tab' or 'T0'
line_val = np.array([21, 22, 23])				# values of each line


Nval = len(line_val)
max_val = {'p':1, 'T0':12, 'Tab':24}
label = {'p':r'$p$', 'T0':r'$T_0$', 'Tab':r'$T_{AB}$'}
val_index = np.array([int(ab_res * val / max_val[line_var]) for val in line_val])
print(val_index)

####################
## IMPORTING DATA ##
####################
lag_opt = np.zeros([Nval, 400])
del_opt = np.zeros([Nval, 400])

# competition optimals
# lag_comp = np.zeros([Nval, 100])
# del_comp = np.zeros([Nval, 100])

i = 0
for index in val_index:
	if line_var == 'p':
		lag_opt[i] = np.loadtxt(data_path + 'high_resolution/optimal_lag-'+ time)[index]
		del_opt[i] = np.loadtxt(data_path + 'high_resolution/optimal_delta-' + time)[index]
		# lag_comp[i] = np.loadtxt(data_path + 'competition_two_species/competition_lag-' + time)[:, index]
		# del_comp[i] = np.loadtxt(data_path + 'competition_two_species/competition_lagdelta-' + time)[:, index]
		color = ['dodgerblue', 'blue', 'black']
	       
	else:
		lag_opt[i] = np.loadtxt(data_path + 'high_resolution/optimal_lag-'+ time)[:,index]
		del_opt[i] = np.loadtxt(data_path + 'high_resolution/optimal_delta-' + time)[:,index]
		# lag_comp[i] = np.loadtxt(data_path + 'competition_two_species/competition_lag-' + time)[index]
		# del_comp[i] = np.loadtxt(data_path + 'competition_two_species/competition_lagdelta-' + time)[index]
		color = ['orange', 'red', 'maroon', 'blue']

	i += 1

x_arr = np.linspace(0, max_val[x_var], 400)
# x_comp = np.linspace(0, max_val[x_var], 100)



####################
## PLOTTING LINES ##
####################
fig, ax = plt.subplots(1, 2, figsize=(12, 5.5), sharex=True)
ax[0].set(xlabel=label[x_var], ylabel=r"$\lambda^*$")
ax[1].set(xlabel=label[x_var], ylabel=r"$\delta^*$", ylim=(0, 0.075))
fig.suptitle(r"Optimal persistence strategy for" + label[time_label] + " = " + str(time_value))

# plotting
for i in range(Nval):
	# ax[0].plot(x_comp, lag_comp[i], 'o', mew=1.5, alpha=1, color=color[i], fillstyle='none',
	#            label=r"$p$ = " + str(p[i % Np]))
	# ax[1].plot(x_comp, del_comp[i], 'o', mew=1.5, alpha=1, color=color[i], fillstyle='none',
	#            label=r"$p$ = " + str(p[i % Np]))

	mask0 = (lag_opt[i] < 1)
	mask1 = (lag_opt[i] > 0.1)
	ax[0].plot(x_arr[mask0], lag_opt[i][mask0], '-', alpha=0.7, lw=4, color=color[i], fillstyle='none',  label=label[line_var] + " = " + str(line_val[i]))
	ax[0].plot(x_arr[mask1], lag_opt[i][mask1], '-', alpha=0.7, lw=4, color=color[i], fillstyle='none')
	# ax[0].plot(x_arr, lag_opt[i], '-', color=color[i], fillstyle='none')

	ax[1].plot(x_arr[mask0], del_opt[i][mask0], '-', alpha=0.7, lw=4, color=color[i], fillstyle='none',  label=label[line_var] + " = " + str(line_val[i]))
	ax[1].plot(x_arr[mask1], del_opt[i][mask1], '-', alpha=0.7, lw=4, color=color[i], fillstyle='none')
	#ax[1].plot(x_arr, del_opt[i], '-', color=color[i], fillstyle='none')


fig.tight_layout(rect=[0, 0, 1, 0.9])
ax[0].legend(loc='upper center', bbox_to_anchor=(1.08, 1.3), ncol=3, fancybox=True, shadow=False)
#ax[1].set(xscale="log", yscale="log", ylim=(10**-4, 1))

# saving
fig.show()
fig.savefig("test.png")
#if saving:
#	fig.savefig("../../../figures/model2/optimal_lines_" + x_var + "-" + time + ".png", dpi=100)
	# fig.savefig("../../../figures/model2/optimal_competitor_lines_" + x_var + "-" + time + ".png", dpi=100)
