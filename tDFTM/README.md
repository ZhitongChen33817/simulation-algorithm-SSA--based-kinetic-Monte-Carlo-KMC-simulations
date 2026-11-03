# Stochastic simulation of sequential assembly on the tDF template (single-chain model)

This folder contains the Gillespie SSA code used for the one-dimensional
template model: oligomer B binds sequentially to the tDF template (A),

```
A + B  <->  AB  <->  AB2  <-> ... <-> AB(n)
```

with `n <= 9`. The network contains 11 species (A, B, AB ... AB9) and 18
unidirectional reaction channels (9 reversible binding steps).

---

## Contents

```
tDFTM/
├── README.md                        this file
├── code/
│   └── tDFTM_SSA_single_chain.py    SSA engine for the single-chain model
└── input/
    └── seq_tDF_9unit.csv            normalised rate constants (one per row)
```

## Model and rate constants

Rate constants are read from `input/seq_tDF_9unit.csv` (column 0). They are
already normalised to molecule-count units for a box of volume `V`:

```
k_on'  = Kon / (N_A * V)     second-order: a = k_on' * n_A * n_B
k_off' = Koff                first-order:  a = k_off' * n_AB
```

`V` only sets the noise amplitude and the number of events per unit time; it
does not change ensemble-averaged concentrations.

Default settings in the script:

| Parameter | Value |
|---|---|
| `n` (maximum number of B units) | 8 |
| `[A]0` (tDF template) | 100 nM |
| `[B]0` | `n * 150 nM` |
| Box volume `V` | 1e-11 L |
| `t_max` | 1000 s |

## Running

Run from this folder (`tDFTM/`), so that the relative input path resolves:

```bash
python code/tDFTM_SSA_single_chain.py
```

The script prints species/reaction counts, the normalised rate constants, and
TDF conservation. To simulate the 9-unit case, set `n = 9` and adjust `con_B`.

## Output

`onehalf_KMC_TDF_100nM_data_8unit.csv`

One row per saved time point, columns:

```
time, A, B, AB, AB2, ..., AB9
```

Values are molecule counts in the simulation box; multiply by `1/(N_A*V)` to
obtain molar concentrations. The sampling interval is adaptive (0.025 s early,
then 1 s, 10 s, 100 s) to resolve the fast initial phase without inflating the
file.

## Reproducibility

Set `random.seed(...)` (or `SEED` if exposed) to fix the trajectory. Ensemble
averages converge as `1/sqrt(N)`; increasing `V` or repeating the run reduces
stochastic scatter without changing the mean behaviour.
