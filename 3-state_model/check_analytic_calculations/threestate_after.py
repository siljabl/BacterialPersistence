import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

sys.path.append("../src")
from differential_equations import S0, f, γ, ode_kill_single
from analytical_calculations import a_b, ap_bp, solve_constants

# setting parameters
d_0 = 1000
λ_d = 27.5099
λ_r =  25.01 
δ = 0.1

t_max = 10
T_0 = 5
T = 5

a, b = a_b(λ_r, δ)
ap, bp = ap_bp(λ_r, δ)
c = 1/λ_d

t_max = T+10


def BAC(a,ap,b,bp,c,T_0, T):
	bac_args = [0, 0, 0, a, b, c, ap, bp]
	ab_args  = [0, T_0, T-T_0]
	
	B, A, C = solve_constants(bac_args, ab_args, stage="post")
	return B, A, C



# functions for computing bacteria density
def d(t, λ_d, λ_r, δ):
	c = 1/λ_d

	return d_0 * np.exp(-c*t)

def g(t, λ_d, λ_r, δ):
	a, b = a_b(λ_r, δ)
	ap, bp = ap_bp(λ_r, δ)
	c = 1/λ_d

	B, A, C = BAC(a,ap,b,bp,c,T_0,T)

	return B*np.exp(b*(t-T)) + A*np.exp(-a*(t-T)) + C*np.exp(-c*t)


def r(t, λ_d, λ_r, δ):
	a, b = a_b(λ_r, δ)
	c = 1/λ_d
	B, A, C = BAC(a,ap,b,bp,c,T_0,T)

	B *=  (1-b)/b
	A *= -(a+1)/a
	C *= (a-b-c-a*b) / (a*b)
	C -= c*d_0 / (a*b)

	return B*np.exp(b*(t-T)) + A*np.exp(-a*(t-T)) + C*np.exp(-c*t)


def n(t, λ_d, λ_r, δ):
	a, b = a_b(λ_r, δ)
	ap, bp = ap_bp(λ_r, δ)
	c = 1/λ_d
	
	B, A, C = BAC(a,ap,b,bp,c, T_0, T)

	return B*np.exp(b*(t-T))/b - A*np.exp(-a*(t-T))/a + (a-b-c)*C*np.exp(-c*t)/(a*b) + (a*b-c)*d_0*np.exp(-c*t)/(a*b)


def ode_grow(t, n, λ_d, λ_r, δ):
    n[n < 0] = 0                                    # avoid negative populations

    dn_dt = np.zeros_like(n)
    dn_dt[0] = -n[0] / λ_d
    dn_dt[1] = n[0] / λ_d + n[2] / λ_r + n[1] * (1-δ)
    dn_dt[2] = -n[2] / λ_r + n[1] * δ

    return dn_dt

# PLOTTING
t = np.linspace(0, t_max, 100)

n0_d = [d(T, λ_d, λ_r, δ), g(T, λ_d, λ_r, δ), r(T, λ_d, λ_r, δ), n(T, λ_d, λ_r, δ)]
sol = solve_ivp(ode_grow, [T, t_max], n0_d, args=(λ_d, λ_r, δ), max_step=0.001)


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
fig.savefig("Check_analytic_population_after_AB.png")

