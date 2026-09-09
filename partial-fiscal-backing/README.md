# Partial Fiscal Backing

Dynare implementations of New Keynesian models with partial fiscal backing.

## Models

- `rank_pfb_mechanism.mod`: a stylized New Keynesian model with a monetary-led economy and a fiscal-led shadow economy.
- `tank_pfb_mechanism.mod`: a TANK New Keynesian model with Ricardian and hand-to-mouth households, capital, long-duration nominal government debt, fiscal policy rules, and a fiscal-led shadow economy.
- `tank_steady_state.m`: auxiliary function used by the TANK model to compute the deterministic steady state.

## Requirements

The models require MATLAB and Dynare. Run each model from its own directory using:

```matlab
dynare rank_pfb_mechanism.mod
```

or:

```matlab
dynare tank_pfb_mechanism.mod
```

The TANK specification automatically produces impulse-response figures and result files.

## References

- Smets, F., and Wouters, R. *Fiscal Backing, Inflation and US Business Cycles*. [CEPR Discussion Paper](https://cepr.org/publications/dp19791).
- Bianchi, F., Faccini, R., and Melosi, L. *A Fiscal Theory of Persistent Inflation*. [NBER Working Paper](https://www.nber.org/system/files/working_papers/w30727/w30727.pdf).

Both codes are part of the dissertation project and applications available at the [University of São Paulo Repository](https://repositorio.usp.br/item/003316180).

## Author

André Chadú
