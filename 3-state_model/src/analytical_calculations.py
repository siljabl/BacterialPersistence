# Coupled ODEs that make up the model
import numpy as np
from differential_equations import f, S0, γ


##########################
## a-b scheme equations ##
##########################
def compute_a_and_b(λ_p, δ):
    '''
    Computing a, b as defined in thesis
    '''
    a, b = 1 / λ_p, np.ones_like(δ)

    D = np.sqrt((λ_p * (1 - δ) - 1) ** 2 + 4 * λ_p)

    a = (D - (λ_p * (1 - δ) - 1)) / (2 * λ_p)
    b = (D + (λ_p * (1 - δ) - 1)) / (2 * λ_p)

    return a, b



def compute_ap_and_bp(λ_p, δ):
    '''
    Computing ap, bp as defined in thesis
    '''
    ap, bp = γ * np.ones_like(δ), 1 / λ_p

    D = np.sqrt((λ_p * (γ + δ) + 1) ** 2 - 4 * γ * λ_p)

    ap =  (D + (λ_p * (γ + δ) + 1)) / (2 * λ_p)
    bp = -(D - (λ_p * (γ + δ) + 1)) / (2 * λ_p)

    return ap, bp



##########################
## Analytical equations ##
##########################
def solve_constants(eq_params, ab_params, stage):
    '''
    Computing constants from equations for growing populations.
    '''
    a  = eq_params['a']
    b  = eq_params['b']
    c  = eq_params['c']
    ap = eq_params['ap']
    bp = eq_params['bp']

    T0  = ab_params['T0']
    Tab = ab_params['Tab']
    T   = ab_params['T']

    ebT0, eaT0, ecT0 = np.exp(b*T0), np.exp(-a*T0), np.exp(-c*T0)
    ebT, eaT, ecT = np.exp(-bp*Tab), np.exp(-ap*Tab), np.exp(-c*Tab)

    B0 =  (c*f*S0 / (a+b)) * (b*(a+1)) / (b+c)
    A0 = -(c*f*S0 / (a+b)) * (a*(1-b)) / (a-c)
    C0 = -(c*f*S0 / (b+c)) * ((c-a*b)) / (c-a)

    #B0[1:] *= (b*(a+1))[1:] / (b[1:]+c[1:])
    #A0[1:] *= (a*(1-b))[1:] / (a[1:]-c[1:])
    #C0[1:] *= ((c-a*b))[1:] / (c[1:]-a[1:])
    if stage == 'pre':
        return B0, A0, C0

    Bp =  ((a-bp)*B0*ebT0 - (b+bp)*A0*eaT0) / (ap-bp)
    Ap = -((a-ap)*B0*ebT0 - (b+ap)*A0*eaT0 + C0*(a-ap)*(b+ap)*ecT0/(ap-c)) / (ap-bp)
    Cp = c*f*S0 / (ap-c)

    #Bp[1:] += (C0*(a-bp)*(b+bp)*ecT0 / (ap-bp))[1:] / (bp[1:]-c[1:])
    #Cp[1:] *= (c-a*b)[1:] / (c[1:]-bp[1:])

    Bp += (C0*(a-bp)*(b+bp)*ecT0 / (ap-bp)) / (bp-c)
    Cp *= (c-a*b) / (c-bp)
    if stage == 'ab':
        return Bp, Ap, Cp
    
    ebT, eaT, ecT = np.exp(-bp*(T-T0)), np.exp(-ap*(T-T0)), np.exp(-c*T)

    B = ((b+ap)*Bp*ebT + (b+bp)*Ap*eaT + Cp*(b+ap)*(b+bp)*ecT/(b+c)) / (a+b)
    A = ((a-ap)*Bp*ebT + (a-bp)*Ap*eaT)/ (a+b)
    C = C0

    #A[1:] += (Cp*(a-ap)*(a-bp)*ecT / (a+b))[1:] / (a[1:]-c[1:])
    A += (Cp*(a-ap)*(a-bp)*ecT / (a+b)) / (a-c)
    if stage == 'post':
        return B, A, C

    else:
        return 0, 0, 0


