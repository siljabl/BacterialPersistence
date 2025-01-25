import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

sys.path.append("../src")
from differential_equations import S0, f, gamma, ode_kill_single, ode_grow_single
from analytical_calculations import compute_a_and_b, compute_ap_and_bp, solve_constants

# setting parameters
d_0 = 1000
λ =  25.01 
δ = 0.1

t_max = 20
T_0 = 5
T = 2

a, b = compute_a_and_b(λ, δ)
ap, bp = compute_ap_and_bp(λ, δ)

t_max = T+10


def BAC(a,ap,b,bp,T_0, T):
	bac_args = {'λ':0, 'δ':0, 'a':a, 'b':b, 'ap':ap, 'bp':bp}
	ab_args  = {'T0':T_0, 'Tab':T-T_0, 'T':0}
	
	B, A = solve_constants(bac_args, ab_args, stage="post")
	return B, A



# functions for computing bacteria density
def d(t, λ, δ):
	a, b = compute_a_and_b(λ, δ)
	ap, bp = compute_ap_and_bp(λ, δ)

	B, A = BAC(a,ap,b,bp,T_0,T)

	return ((1-b)/b)*B*np.exp(b*(t-T)) - ((a+1)/a)*A*np.exp(-a*(t-T))


def g(t, λ, δ):
	a, b = compute_a_and_b(λ, δ)
	ap, bp = compute_ap_and_bp(λ, δ)

	B, A = BAC(a,ap,b,bp,T_0,T)

	return B*np.exp(b*(t-T)) + A*np.exp(-a*(t-T))


# PLOTTING
t = np.linspace(0, t_max, 100)


n0_d = [d(T, λ, δ), g(T, λ, δ)]
sol = solve_ivp(ode_grow_single, [T, t_max], n0_d, args=(λ, δ), max_step=0.001)


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
	ax[1,i].set(yscale="linear")
	ax[0,i].legend()
fig.savefig("Check_analytic_population_after_AB.png")

