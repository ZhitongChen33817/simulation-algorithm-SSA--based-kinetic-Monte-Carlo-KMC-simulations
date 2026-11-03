# Stochastic simulation code for the tDFTM two-layer reaction network

This package contains the Python implementation of the Gillespie stochastic
simulation algorithm (SSA) used to model the sequential assembly of oligomer B
onto the tDF template (tDFTM) and the secondary binding of oligomer C to the
C-sites of the growing B chain. It reproduces the parameter scans and figures
reported in the manuscript.

---

## 1. Model

The tDF template provides `max_B` primary B-binding sites. B assembles
sequentially:

```
TDF + B  <->  TDF-B1  <-> ... <-> TDF-B(max_B)
```

Each bound B creates a secondary C-site, but only once the main chain has
advanced far enough: a C oligomer bound at position `j` requires `j <= i - 2`,
where `i` is the number of bound B units. This constraint generates a
two-dimensional species space

```
TDF-B(i)-C(j),   0 <= i <= max_B,   0 <= j <= max(0, i-2)
```

Every species participates in two reversible reactions (one B addition/removal
and, when allowed, one C addition/removal), so the network has
`2 * N_species` reaction channels. Species names, index maps and reaction
stoichiometries are generated programmatically from these rules
(`build_species`, `build_reactions` in `code/tDFTM_network_SSA.py`).

### Simulation algorithm

Gillespie's direct method is used. For every step:

1. propensities `a_k = k_on * n_X * n_Y` (bimolecular) or `a_k = k_off * n_X`
   (unimolecular) are computed for all channels;
2. `a0 = sum(a_k)`, `tau = -ln(u1) / a0`;
3. a channel is selected with `sum(a_j) >= u2 * a0`.

Rate constants follow the standard concentration-to-molecule-count
normalisation for a box of volume `V`:

```
k_on'  = Kon / (N_A * V)      [s^-1 per molecule^-1]
k_off' = Koff                  [s^-1]
Kd     = Koff / Kon            [M]
```

`V` sets the noise amplitude and the number of events per unit time; it does
not change ensemble-averaged concentrations. The manuscript data were produced
with `V = 1e-11 L`.

### Parameters used for the manuscript

| Parameter | Value |
|---|---|
| Avogadro constant `N_A` | 6.022e23 mol^-1 |
| Box volume `V` (production runs) | 1e-11 L |
| `[TDF]0` | 100 nM |
| `max_B` | 3 |
| `[B]0` per site | 1 - 1000 nM (= ratio x 100 nM) |
| `ratio` ([B]/max_B : TDF) | 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 5, 7.5, 10 |
| `Kon_B` | 1.66e6 M^-1 s^-1 |
| `Koff_B` | 3.32e-5 s^-1 (`Kd_B` = 20 pM) |
| `Kd_C` | 50 nM, 100 nM, 200 nM, 500 nM, 1 uM, 2 uM, 5 uM, 10 uM, 20 uM, 50 uM |
| `Kon_C` | 5.0e2 M^-1 s^-1 |
| `Koff_C` | `Kd_C * Kon_C` |
| Simulation time `t_max` | 10000 s |

For each `(ratio, Kd_C)` pair the network is simulated once with `max_B = 3`;
the time series is written in full so that kinetics and endpoint distributions
can both be analysed.

---

## 2. Package contents

```
EM/
├── README.md                  this file
├── requirements.txt           Python dependencies
└── code/
    ├── tDFTM_network_SSA.py   core model + Gillespie engine + single run
    └── scan.py                batch scan over ratio x Kd_C (manuscript data)
```

All rate constants are defined as constants at the top of
`code/tDFTM_network_SSA.py`; no external parameter file is required.

---

## 3. Requirements

Python 3.8 or newer, with:

```bash
pip install -r requirements.txt
```

Dependencies: `numpy`, `pandas`, `matplotlib` (Agg backend, no display needed).
`xlsxwriter` is optional and only needed if the Excel-export helper scripts are
added.

---

## 4. Running the code

### 4.1 Single simulation

```bash
cd code
python tDFTM_network_SSA.py        # defaults: max_B = 9, t_max = 1000 s
python tDFTM_network_SSA.py 5000   # override t_max (s)
```

Edit the constants at the top of `tDFTM_network_SSA.py` (`max_B`, `con_B`,
`Kd_C`, `V`, `SEED`) to change the condition.

### 4.2 Full parameter scan (manuscript data set)

```bash
cd code
python scan.py            # uses t_max = 10000 s
python scan.py 1000       # optional: override t_max
```

This loops over `ratio_list` (outer), `max_B_list`, and `Kd_C_list` (inner) and
writes one directory (`source_data_1030/`) containing:

| File | Content |
|---|---|
| `kinetics_TDF1_B{ratio}_maxB{n}_KdC{k}nM.csv` | full time series of every species (one row per saved time point) |
| `TDF1_B{ratio}_scan_summary.csv` | one row per `Kd_C`: `avg_primary_B`, `avg_C_per_TDF`, `frac_with_C`, `frac_bare_TDFBmax`, `free_B`, `n_events`, `t_end` |
| `TDF1_B{ratio}_scan_final_states.csv` | long format endpoint distribution of all species |

### 4.3 Output format

`kinetics_*.csv` columns:

```
time, B, TDF, TDFB1, TDFB2, ..., TDFBn, TDFBnC1, ...
```

All species columns are **molecule counts** in the simulation box (multiply by
`1/(N_A*V)` for molar concentrations). The sampling interval is adaptive:
0.025 s below 30 s, 1 s until 200 s, 10 s until 2000 s, then 100 s, which keeps
the early fast phase resolved while limiting file size.

### 4.4 Figures and tables

Figures in the manuscript were generated from the `kinetics_*.csv` files written
by `scan.py`. Each file contains the full species-resolved time series, so both
the endpoint distribution (last row) and the kinetic curves (any row/column
pair) can be reproduced directly. The plotting routines are standard
matplotlib scripts and are not included in this archive.

---

## 5. Reproducibility notes

- `SEED = None` in `scan.py` / `tDFTM_network_SSA.py` gives a fresh random
  stream per run. Set an integer to make a single trajectory reproducible.
- Ensemble averages converge as `1/sqrt(N)`; increase `V` (or the number of
  independent repeats) to reduce stochastic scatter without changing the
  mean behaviour.
- The C-site rule `j <= i - 2` is enforced when the species/reaction tables are
  built, so no reaction can create a forbidden state.
- Species counts are stored as integers; downstream scripts convert them to
  fractions of total TDF (`sum of all TDF-containing species`).

---

## 6. Manuscript-ready Methods paragraph (draft)

> **Stochastic simulation of the two-layer assembly network.** The reaction
> network was simulated with Gillespie's direct stochastic simulation
> algorithm (SSA). The tDF template was modelled as a chain of `max_B` primary
> B-binding sites that are occupied sequentially; each bound B unit creates a
> secondary C-site, with the geometric constraint that a C oligomer can occupy
> position `j` only when at least `i >= j + 2` B units are already bound. Every
> species `TDF-B(i)-C(j)` therefore participates in a reversible B
> addition/removal step and, where allowed, a reversible C addition/removal
> step. Bimolecular propensities were computed as `Kon/(N_A V) * n_X * n_Y`
> and unimolecular propensities as `Koff * n_X`, with `N_A` the Avogadro
> constant and `V = 1e-11 L` the simulation volume. The B-site kinetics used
> `Kon_B = 1.66e6 M^-1 s^-1` and `Koff_B = 3.32e-5 s^-1` (`Kd_B` = 20 pM); the
> C-site affinity was varied over `Kd_C` = 50 nM - 50 uM with
> `Koff_C = Kd_C * Kon_C`. Initial concentrations were `[TDF]0` = 100 nM and
> `[B]0` = `ratio * 100 nM` per site, with `ratio` spanning 0.01-10. Each
> condition was propagated to `t_max` = 10,000 s, and the full species-resolved
> time series was recorded. All simulations were performed in Python 3 using
> custom code (provided as Supplementary Code).

---

## 7. Data and code availability

The simulation code is included in this package. The complete scan output
(`source_data_1030/`, 190 files: 19 ratios x 10 `Kd_C` values, plus per-ratio
summary tables) is available from the corresponding author upon request, or can
be regenerated with `python scan.py 10000`.
