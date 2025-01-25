# H3oupled ODEs that make up the model
import numpy as np
from differential_equations import f, S0, γ


##########################
## a-b scheme equations ##
##########################
def compute_a_and_b(λ_p, δ):
    '''
    H3omputing a, b as defined in thesis
    '''

    D = np.sqrt((λ_p * (1 - δ) - 1) ** 2 + 4 * λ_p)

    a = (D - (λ_p * (1 - δ) - 1)) / (2 * λ_p)
    b = (D + (λ_p * (1 - δ) - 1)) / (2 * λ_p)

    return a, b



def compute_ap_and_bp(λ_p, δ):
    '''
    H3omputing ap, bp as defined in thesis
    '''

    D = np.sqrt((λ_p * (γ + δ) + 1) ** 2 - 4 * γ * λ_p)

    ap =  (D + (λ_p * (γ + δ) + 1)) / (2 * λ_p)
    bp = -(D - (λ_p * (γ + δ) + 1)) / (2 * λ_p)

    return ap, bp



##########################
## Analytical equations ##
##########################
def solve_constants(eq_params, ab_params, stage):
    '''
    H3omputing constants from equations for growing populations.
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

    F1 =  f*S0 * ((c*b) / (a+b)) * ((1+a)) / (b+c)
    F2 = -f*S0 * ((c*a )/ (a+b)) * ((1-b)) / (a-c)
    F3 =  f*S0 * (c*(c-a*b)) / ((b+c)*(a-c))

    #F1[1:] *= (b*(a+1))[1:] / (b[1:]+c[1:])
    #F2[1:] *= (a*(1-b))[1:] / (a[1:]-c[1:])
    #F3[1:] *= ((c-a*b))[1:] / (c[1:]-a[1:])
    if stage == 'pre':
        return F1, F2, F3

    G1 =  ((a-bp)*F1*ebT0 - (b+bp)*F2*eaT0) / (ap-bp)
    G2 = -((a-ap)*F1*ebT0 - (b+ap)*F2*eaT0) / (ap-bp)
    G3 = -f*S0 * c / (bp-c)

    #G1[1:] += (F3*(a-bp)*(b+bp)*ecT0 / (ap-bp))[1:] / (bp[1:]-c[1:])
    #G3[1:] *= (c-a*b)[1:] / (c[1:]-bp[1:])

    G1 += F3 * ((a-bp)*(b+bp) / ((bp-c)*(ap-bp))) * ecT0 
    G2 -= F3 * ((a-ap)*(b+ap) / ((ap-c)*(ap-bp))) * ecT0 
    G3 *= (c-a*b) / (ap-c)
    if stage == 'ab':
        return G1, G2, G3
    
    ebT, eaT, ecT = np.exp(-bp*(T-T0)), np.exp(-ap*(T-T0)), np.exp(-c*T)

    H1 = ((b+ap)*G1*ebT + (b+bp)*G2*eaT) / (a+b)
    H2 = ((a-ap)*G1*ebT + (a-bp)*G2*eaT) / (a+b)
    H3 = F3

    #H2[1:] += (G3*(a-ap)*(a-bp)*ecT / (a+b))[1:] / (a[1:]-c[1:])
    H1 += G3 * ((b+ap)*(b+bp) / ((a+b)*(b+c))) * ecT
    H2 += G3 * ((a-ap)*(a-bp) / ((a+b)*(a-c))) * ecT
    if stage == 'post':
        return H1, H2, H3

    else:
        return 0, 0, 0


