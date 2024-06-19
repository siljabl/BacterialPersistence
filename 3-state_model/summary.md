# Summary on 3 state simulations
## "Method"
I originally computed the optimal parameters as in the two state case, that is by computing the averate time it would take to consume all the nutrients ($T_S$). However, to isolate $T_S$ I had to make an assumption that I on closer inspection did not trust. Now I have redone the calculations, and compute the optimal parameters without the previous assumption, but instead by numerically determining $T_S$. For a given set of antibiotic parameters (p, $T_0$, $T_{ab}$) I determine $T_S$ for every set of bacterial parameters ($\lambda_d$, $\lambda_r$, $\delta$). The optimal combination of ($\lambda_d$, $\lambda_r$, $\delta$) is the one that minimizes $T_S$.

The competition average parameters are computed by evolving several species according to the differential equations and using a solver to find $T_S$ for 10 000 consecutive cycles.

The mutation simulations are done like the competition simulations, but with a mutation rate between the different species. Every simulation is started from a single species with a specific set of bacterial parameters (min and max?)

- What about mutation from $\lambda_d$=0? Create exception?
- Extinction?

## Coupled nutrients and antibiotics, $T_0 = 0$
![Optimal parameters for $T_0=0$](figs/single_optimal/optimal_heatmap_T0_0.png)

I still get same result as before: optimal strategy is either only triggered persistence, or only spontaneous persistence. The result is confirmed by a competition simulation, with competition between many species with different combinations of parameters.
Having very small fraction of spontaneous persistence affects the growth very little, however after very many cycles the loss is significant.

![Interspecies competition at $T_0=0$, $T_{ab}=12$](figs/competition_average/average_parameters-T0_0-T_12.png)

[comment]: <> (![The bacterial parameters from the competition simulation for $T_0=0$]())

![10 example feast-famine cycles]()


## Decoupled nutrients and antibiotics, $T_0 > 0$
![Optimal parameters for $T_{AB}=12$](figs/single_optimal/optimal_heatmap_Tab_12.png)
Also when the antibiotics are decoupled from the addition of nutrients the two strategies are separated.

![Interspecies competition at $T_0=2$, $T_{ab}=12$]()

![10 example feast-famine cycles]()


## Mutation
Redo competition simulations with mutation.

### To do
1) rerun heatmaps high resolution
2) rerun competition to confirm
3) Modify to add mutation
4) Run mutation



