import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

sys.path.append("../src")
from differential_equations import S0, f, ode_grow_single
from analytical_calculations import compute_a_and_b, compute_ap_and_bp, solve_constants

# setting parameters
d_0 = f*S0
λ_d = 2
λ_r = 10
δ = 0.1

a, b = compute_a_and_b(λ_r, δ)
c = 1/λ_d


def BAC(a,b,c):
	bac_args = {'λ_d':0, 'λ_d':0, 'δ':0, 'a':a, 'b':b, 'c':c, 'ap':0, 'bp':0}
	ab_args  = {'T0':0, 'Tab':0, 'T':0}
	
	B, A, C = solve_constants(bac_args, ab_args, stage="pre")
	return B, A, C


# functions for computing bacteria density
def d(t, λ_d, λ_r, δ):
	a, b = compute_a_and_b(λ_r, δ)
	c = 1/λ_d

	return d_0 * np.exp(-c*t)

def g(t, λ_d, λ_r, δ):
	a, b = compute_a_and_b(λ_r, δ)
	c = 1/λ_d
	B, A, C = BAC(a, b, c)

	return B*np.exp(b*t) + A*np.exp(-a*t) + C*np.exp(-c*t)


def r(t, λ_d, λ_r, δ):
	a, b = compute_a_and_b(λ_r, δ)
	c = 1/λ_d
	B, A, C = BAC(a, b, c)

	B *=  (1-b)/b
	A *= -(a+1)/a
	C *= -(a+1)*(1-b) / (c-a*b)
	return B*np.exp(b*t) + A*np.exp(-a*t) + C*np.exp(-c*t)


# PLOTTING
t_max = 50
t = np.linspace(0, t_max, 100)

n0_d = [d_0, 0, 0]
sol = solve_ivp(ode_grow_single, [0, t_max], n0_d, args=(λ_d, λ_r, δ), max_step=0.001)

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
ax[0,2].plot(sol.t, r(sol.t, λ_d, λ_r, δ), label="r(t)")
ax[0,2].plot(sol.t, sol.y[2], "--", label="ode")
ax[1,2].plot(sol.t, sol.y[2]-r(sol.t, λ_d, λ_r, δ))

for i in range(3):
	ax[0,i].set(yscale="linear")
	ax[0,i].legend()

fig.savefig("Check_analytic_population_before_AB.png")

def analytical_fitness(bac_args, ab_args):
    _, _, a, b, ap, bp = bac_args
    p, T0, Tab = ab_args

    # consumption time without antibiotics
    Ts = (1 / b) * np.log((a + b) / (a * b * f))    

    prefactor = a * b * f / ((ap - bp) * (a + b) ** 2)
    compute_a_and_bp = a - bp
    a_ap = a - ap
    b_ap = b + ap
    b_bp = b + bp

    exp_aT0 = np.exp(-a * T0)
    exp_bT0 = np.exp(b * T0)
    exp_apT = np.exp(-ap * Tab)
    exp_bpT = np.exp(-bp * Tab)

    gT0 = a * b * f * (exp_bT0 - exp_aT0) / (a + b)
    gT = prefactor * (a + b) * (compute_a_and_bp * exp_bT0 + b_bp * exp_aT0) * exp_bpT - \
         prefactor * (a + b) * (a_ap * exp_bT0 + b_ap * exp_aT0) * exp_apT

    b_term = b_ap * (compute_a_and_bp * exp_bT0 + b_bp * exp_aT0) * exp_bpT - \
             b_bp * (a_ap * exp_bT0 + b_ap * exp_aT0) * exp_apT

    # consumption time with antibiotics
    Ts_ab = (1 / b) * np.log((1 + gT - gT0) / (prefactor * b_term))
    
    return Ts_ab

# setting parameters
d_0 = f*S0
λ_d = 0.01
λ_r = 0.01
δ = 0
T_0=0
T=24

a, b = compute_a_and_b(λ_r, δ)
ap, bp = compute_ap_and_bp(λ_r, δ)
c = 1/λ_d

bac_args = [λ_d, δ, a, b, ap, bp]
ab_args = [0, T_0, T]

print(analytical_fitness(bac_args, ab_args))