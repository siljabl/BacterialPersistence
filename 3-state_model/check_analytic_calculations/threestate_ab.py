import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

sys.path.append("../src")
from differential_equations import S0, f, γ, ode_kill_single
from analytical_calculations import a_b, ap_bp, solve_constants

# setting parameters
d_0 = 1000
λ_d = 2
λ_r = 1.5
δ = 0.1

t_max = 10
T_0 = 0

a, b = a_b(λ_r, δ)
ap, bp = ap_bp(λ_r, δ)
c = 1/λ_d


def BAC_p(a,ap,b,bp,c,T_0):
	bac_args = [0, 0, 0, a, b, c, ap, bp]
	ab_args  = [0, T_0, 0]
	
	Bp, Ap, Cp = solve_constants(bac_args, ab_args, stage="ab")
	return Bp, Ap, Cp



# functions for computing bacteria density
def d(t, λ_d, λ_r, δ):
	c = 1/λ_d

	return d_0 * np.exp(-c*t)


def g(t, λ_d, λ_r, δ):
	a, b = a_b(λ_r, δ)
	ap, bp = ap_bp(λ_r, δ)
	c = 1/λ_d

	Bp, Ap, Cp = BAC_p(a,ap,b,bp,c,T_0)

	return Bp*np.exp(-bp*(t-T_0)) + Ap*np.exp(-ap*(t-T_0)) + Cp*np.exp(-c*t)


def r(t, λ_d, λ_r, δ):
	a, b = a_b(λ_r, δ)
	ap, bp = ap_bp(λ_r, δ)
	c = 1/λ_d

	Bp, Ap, Cp = BAC_p(a,ap,b,bp,c,T_0)

	Bp *= (ap - a*b) / (a*b)
	Ap *= (bp - a*b) / (a*b)
	Cp *= (ap+bp-a*b-c) / (a*b)
	Cp -= c*d_0 / (a*b)

	return Bp*np.exp(-bp*(t-T_0)) + Ap*np.exp(-ap*(t-T_0)) + Cp*np.exp(-c*t)


# PLOTTING
t = np.linspace(0, t_max, 100)

n0_d = [d(T_0, λ_d, λ_r, δ), g(T_0, λ_d, λ_r, δ), r(T_0, λ_d, λ_r, δ)]
sol = solve_ivp(ode_kill_single, [T_0, t_max], n0_d, args=(λ_d, λ_r, δ), max_step=0.001)

fig, ax = plt.subplots(2, 3, figsize = (12,7))
ax[0,0].set(ylabel="Population")
ax[1,0].set(ylabel="Error")

ax[0,0].set(title="Triggered dormant")
ax[0,0].plot(sol.t, d(sol.t, λ_d, λ_r, δ), label="d(t)")
ax[0,0].plot(sol.t, sol.y[0], "--", label="ode")
ax[1,0].plot(sol.t, sol.y[0]-d(sol.t, λ_d, λ_r, δ))

ax[0,1].set(title="Growing")
ax[0,1].plot(sol.t, g(sol.t, λ_d, λ_r, δ), label="g(t)")
ax[0,1].plot(sol.t, sol.y[1], "--", label="ode")
ax[1,1].plot(sol.t, sol.y[1]-g(sol.t, λ_d, λ_r, δ))

ax[0,2].set(title="Spontaneous persistent")
ax[0,2].plot(sol.t, r(sol.t, λ_d, λ_r, δ), label="p(t)")
ax[0,2].plot(sol.t, sol.y[2], "--", label="ode")
ax[1,2].plot(sol.t, sol.y[2]-r(sol.t, λ_d, λ_r, δ))

for i in range(3):
	ax[0,i].set(yscale="linear")
	ax[0,i].legend()

fig.savefig("Check_analytic_population_during_AB.png.png")
