README - Stochastic simulation code for the DNA nanomachine assembly networks
==============================================================================

This archive contains the Gillespie stochastic simulation (SSA / kinetic Monte
Carlo) code used in the manuscript. The code is organised in three
self-contained folders; each folder has its own README.md describing the model,
the input parameters, how to run it and the output format.

Folder 1: 2WJ/
------------------------------------------------------------------
Two-way junction (2WJ) assembly network.
Two orthogonal hybridisation directions are assembled on the same scaffold,
giving a two-dimensional species grid A(i)-B(j) with up to 137 reversible
reactions (274 rate constants). The network can be truncated to a fragment of
interest (e.g. A2B3) by keeping only the relevant rate constants.

  code/    Python_Code_for_SSA_2WJ.py   SSA engine
  input/   seq_ssDNA_9unit.csv          normalised rate constants

Folder 2: tDFTM/
------------------------------------------------------------------
Sequential assembly on the tDF template (single-chain model).
Oligomer B binds sequentially to the tDF template, A + B <-> AB <-> AB2 <->
... <-> AB(n), with n <= 9 (11 species, 18 unidirectional channels).

  code/    tDFTM_SSA_single_chain.py    SSA engine
  input/   seq_tDF_9unit.csv            normalised rate constants

Folder 3: EM/
------------------------------------------------------------------
Two-layer network: tDFTM with a secondary C layer (the data set used for the
manuscript figures). The primary B chain occupies up to max_B sites; a C
oligomer can occupy position j only when at least j+2 B units are bound
(j <= i-2), giving the species space TDF-B(i)-C(j). The folder contains the
scan driver (19 ratios x 10 Kd_C values, max_B = 3, t_max = 10000 s) and the
core engine; all rate constants are set as constants at the top of the engine
file, so no external parameter file is required.

  code/    tDFTM_network_SSA.py         core model + Gillespie engine
  code/    scan.py                      batch scan over ratio x Kd_C
  README.md                             model, parameters, outputs, figure mapping
  requirements.txt

Common requirements
------------------------------------------------------------------
Python 3.8 or newer:

  pip install numpy pandas matplotlib xlsxwriter

(matplotlib is used with the Agg backend, so no display is needed;
xlsxwriter is only used by the Excel export in folder 3 if those helper
scripts are added.)

How to run
------------------------------------------------------------------
Each folder is self-contained. Change into the folder and run its main script:

  cd 2WJ    && python code/Python_Code_for_SSA_2WJ.py
  cd tDFTM  && python code/tDFTM_SSA_single_chain.py
  cd EM     && python code/scan.py 10000

All scripts use Gillespie's direct method with propensities
a = k_on' * n_X * n_Y (bimolecular) or a = k_off' * n_X (unimolecular), with
rate constants normalised for a simulation box of volume V:

  k_on'  = Kon / (N_A * V)
  k_off' = Koff
  Kd     = Koff / Kon

V only sets the noise amplitude and the number of events per unit time; it does
not change ensemble-averaged concentrations.

Output format
------------------------------------------------------------------
Each script writes a CSV with one row per saved time point: the first column is
time (s), subsequent columns are molecule counts of the species in the order
defined by the network. Multiply counts by 1/(N_A*V) for molar concentrations.
The sampling interval is adaptive, so the fast initial phase is resolved while
file size stays manageable.

Reproducibility and data availability
------------------------------------------------------------------
Set the random seed in the script for a reproducible single trajectory.
Ensemble averages converge as 1/sqrt(N); increase V or average independent
repeats to reduce stochastic scatter. The complete scan output of folder 3
(190 kinetic time series, 19 ratios x 10 Kd_C values) can be regenerated with
"python scan.py 10000" and is available from the corresponding author upon
request.
