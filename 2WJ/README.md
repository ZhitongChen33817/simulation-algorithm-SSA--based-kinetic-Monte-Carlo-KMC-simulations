# Stochastic simulation of the two-way junction (2WJ) assembly network

This folder contains the Gillespie SSA-KMC code for the two-way junction
reaction network. Two orthogonal hybridisation directions are assembled on the
same scaffold, giving a two-dimensional species grid

```
A(i)-B(j)
```

with up to 137 reversible reactions (274 rate constants). The full network can
be truncated to a fragment of interest (for example A2B3) by keeping only the
rate constants that belong to that fragment and discarding higher-order
reactions.

---

## Contents

```
2WJ/
├── README.md                        this file
├── code/
│   └── Python_Code_for_SSA_2WJ.py    SSA engine for the 2WJ network
└── input/
    └── seq_ssDNA_9unit.csv          normalised rate constants (274 entries)
```

## Input format

`input/seq_ssDNA_9unit.csv` holds the normalised forward/reverse rate constants
in column 0, one per row (forward and reverse alternating), and the corresponding
species labels in the higher columns (e.g. `A1B1`, `A1B2`, `A2B1`, `A2B2`,
`A2B3`, ...). Each row therefore defines one unidirectional reaction of the
2WJ ladder.

The first rows also record the box volume, `Kon` and `Kd`, which are used only
for reporting:

```
volume (L), Kon (M^-1 s^-1), Kd (M)
```

Rate constants are already normalised to molecule-count units:

```
k_on'  = Kon / (N_A * V)     second-order: a = k_on' * n_X * n_Y
k_off' = Koff                first-order:  a = k_off' * n_X
```

## Running

Run from this folder (`2WJ/`):

```bash
python code/Python_Code_for_SSA_2WJ.py
```

The script reads the CSV above, builds the reaction list, runs Gillespie's
direct method to the requested `t_max`, and writes the time evolution of all
species to a CSV file.

## Output

One row per saved time point; the first column is time (s) and the following
columns are molecule counts of the species in the order defined by the network
table. Convert to molar concentration with `1/(N_A*V)`.

## Notes

- Only reactions with non-zero rate constants in the CSV are active; a subset
  file (e.g. A2B3 only) reproduces the truncated-fragment simulation.
- `V` changes the noise amplitude and the event count per unit time, not the
  ensemble-averaged result.
- Fix the random seed for a reproducible single trajectory; increase `V` or
  average independent repeats to reduce stochastic scatter.
