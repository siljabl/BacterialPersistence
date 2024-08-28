# Summary on 3 state simulations
## "Method"
I originally computed the optimal parameters as in the two state case, that is by analytically computing the average time it would take to consume all the nutrients ($T_S$).
However, to isolate $T_S$ I had to make an assumption that I did not trust on close inspection.
Now I have redone the calculations, and computed the optimal parameters without the previous assumption, but instead by numerically determining $T_S$.
For a given set of antibiotic parameters (p, $T_0$, $T_{ab}$) I determine $T_S$ for every set of bacterial parameters ($\lambda_d$, $\lambda_r$, $\delta$).
The optimal combination of ($\lambda_d$, $\lambda_r$, $\delta$) is the one that minimizes $T_S$.

In addition to the theoretical optimal parameters, I have also computed the competition average parameters.
This is done by evolving several species according to the differential equations and using a solver to find $T_S$ for 20 000 consecutive cycles. The different species have parameters $\lambda_{d/r} \in [0.01, T]$ with $d\lambda = 1$ and $\delta \in [0, 0.05]$ with $d\delta = 0.001$.

The mutation simulations are done like the competition simulations, but with a mutation rate between the different species. Every simulation is started from a single species with a specific set of bacterial parameters (min and max?)

- What about mutation from $\lambda_d$=0? Create exception?
- Extinction?

## Coupled nutrients and antibiotics, $T_0 = 0$
![Optimal parameters for $T_0=0$](figs/single_optimal/optimal_heatmap_T0_0.png)

I still get same result as before: optimal strategy is either only triggered persistence, or only spontaneous persistence (see Fig. 1 and Fig. 5). 
$\lambda_d^*$ is the same as earlyer, i.e. $\lambda_d^* \approx pT$, whereas $\lambda_r^* \approx 0.85T$.
The value of $\delta^*$ is mainly determined by $p$.
I am a bit surprised that $\lambda_r$ is not smaller, since bacteria is entering spontaneous persistence both during and after the antibiotics.
I assume this is balanced by a low $\delta^*$ (though still not as low as experimentaly observed).

![Interspecies competition at $T_0=0$, $T_{ab}=12$](figs/competition_average/average_parameters-T0_0-T_12.png)

The result is confirmed by a competition simulation in Fig. 2, where the dashed lines represent the theoretical optimals from Fig. 1.
For $p > 0.1$ the optimal is to have only triggered persistence, whereas for $p=0.1$ spontaneous persistence is the optimal.
$p=0.3$ is very close to the phase transition, and is therefore fluctuating slightly between the two optimals.

For $\lambda_d$ the competition average is not perfectly consistent with the theoretical optimal, which I think is because the resolution of the parameters in Fig. 1 is much higher (The competition average is much more computationally heavy to compute).

The behaviour of $\delta$ for $p=0.3$ is a bit weird.
What I think happens is that when this weakly bistable system jumps from a low risk state (only spontaneous persistence) to a high risk state (only triggered persistence), it also benefits from the marginal additional protection from having $\delta = \delta_{max} = 0.05$.
With time $\delta$ decreases toward 0, but since $\lambda_r = 0.01$, the penalty for having non-zero $\delta$ is very small, hence the decrease is very slow.
The parameter combination of $\lambda_r \approx 0$ and $\delta > 0$ is probably not very realistic.

The last odd feature of the plot is for $p=0.1$.
Whereas $\lambda_r$ and $\delta$ fluctuate a lot, $\lambda_d$ is not.
For $p=0.1$, antibiotics are so rare that for long periods there are no cycles with antibiotics.
During these periods $\lambda_r\to\lambda_{min}$, but as soon as there is a round of antibiotcs $\lambda_r$ jumps back to the theoretical optimal.
It is not really clear to me why $\delta$ should be increasing in the absence of antibiotics.



## Decoupled nutrients and antibiotics, $T_0 > 0$
![Optimal parameters for $T_{AB}=12$](figs/single_optimal/optimal_heatmap_Tab_12.png)

Also when the antibiotics are decoupled from the addition of nutrients the two strategies are separated (see Fig. 2 and Fig. 6).
Again, $\lambda_d^* \approx pT$, $\lambda_r^*  \approx 0.85T$, and the value of $\delta^*$ is mainly determined by $p$.
I've probably set the upper limit on $\delta$ too low.


![Interspecies competition at $T_0=2$, $T_{ab}=12$](figs/competition_average/average_parameters-T0_2-T_14.png)

I have also run a competition simulation in Fig. 4.
The figure is a bit messy, but still in agreement with the theoretical optimals.
For $p<0.7$ spontaneous persistence is the optimal strategy, and for $p=>0.7$ triggered persistence is the optimal.
However, both $p=0.5$ and $p=0.7$ are close to the phase boundary, with strong fluctuations.
$p=0.7$ shows similar behavior as $p=0.3$ in Fig. 2, but with the decay of $\delta$ being much faster.
I think that is because $\delta^*_{p=0.7}\approx\delta_{max}$ (and I've probably set the upper limit on $\delta$ too low.)
The spikes in $\delta$ where the decay back to the optimal value happens immediately represent flucuations that are not large enough to the system to switch to spontaneous persistence.

$p=0.1, 0.3, 0.5$, behave like $p=0.1$ in Fig. 2, i.e. with fluctuations away from the optimal strategy, but they never switch to triggered persistence.
For $p=0.5$ the optimal stratefgy of triggered persistence has a finite $\lambda_d$, whereas for $p=0.1$ and $p=0.3$ it is $\lambda_{min}\approx0$.
Lastly, I think the fluctuations at $p=0.3$ are smaller than for both $p=0.1$ and 0.5 because the penalty for switchin phase is the highest at $p=0.3$.

## Mutation
In progess


## Rescaled heatmaps
![Optimal parameters for $T_0=0$](figs/single_optimal/optimal_heatmap_T0_0_rescaled.png)

![Optimal parameters for $T_{AB}=12$](figs/single_optimal/optimal_heatmap_Tab_12_rescaled.png)
