import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

sys.path.append("../src")
from differential_equations import S0, f, gamma, ode_kill_single
from analytical_calculations import compute_a_and_b, compute_ap_and_bp, solve_constants

# setting parameters
d_0 = f*S0
λ = 2
δ = 0.1

t_max = 10
T_0 = 0

a,  b  = compute_a_and_b(λ, δ)
ap, bp = compute_ap_and_bp(λ, δ)


def BAC_p(a,ap,b,bp,T_0):
	bac_args = {'λ_d':0, 'δ':0, 'a':a, 'b':b, 'ap':ap, 'bp':bp}
	ab_args  = {'T0':T_0, 'Tab':0, 'T':0}
	
	Bp, Ap = solve_constants(bac_args, ab_args, stage="ab")
	return Bp, Ap



# functions for computing bacteria density
def d(t, λ, δ):
	ap, bp = compute_ap_and_bp(λ, δ)
	a, b = compute_a_and_b(λ, δ)
	Bp, Ap = BAC_p(a,ap,b,bp,T_0)

	return ((ap-a*b)*Bp*np.exp(-bp*(t-T_0)) + (bp-a*b)*Ap*np.exp(-ap*(t-T_0))) / (a*b)


def g(t, λ, δ):
	a,  b  = compute_a_and_b(λ, δ)
	ap, bp = compute_ap_and_bp(λ, δ)

	Bp, Ap = BAC_p(a,ap,b,bp,T_0)

	return Bp*np.exp(-bp*(t-T_0)) + Ap*np.exp(-ap*(t-T_0))



# PLOTTING
t = np.linspace(0, t_max, 100)

n0_d = [d(T_0, λ, δ), g(T_0, λ, δ)]
sol = solve_ivp(ode_kill_single, [T_0, t_max], n0_d, args=(λ, δ), max_step=0.001)

fig, ax = plt.subplots(2, 2, figsize = (8,7))
ax[0,0].set(ylabel="Population")
ax[1,0].set(ylabel="Error")

ax[0,0].set(title="Triggered dormant")
ax[0,0].plot(sol.t, d(sol.t, λ, δ), label="d(t)")
ax[0,0].plot(sol.t, sol.y[0], "--", label="ode")
ax[1,0].plot(sol.t, sol.y[0]-d(sol.t, λ, δ))

ax[0,1].set(title="Growing")
ax[0,1].plot(sol.t, g(sol.t, λ, δ), label="g(t)")
ax[0,1].plot(sol.t, sol.y[1], "--", label="ode")
ax[1,1].plot(sol.t, sol.y[1]-g(sol.t, λ, δ))

for i in range(2):
	ax[0,i].set(yscale="linear")
	ax[0,i].legend()

fig.savefig("Check_analytic_population_during_AB.png.png")
