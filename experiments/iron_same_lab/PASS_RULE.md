# Pass rule: iron(III) rates against the laboratory that built the model, with its one fitted constant held out

Written 2026-10-31, before the iron code existed. De Laat & Gallard (1999) state that their model, with k6 = 2.7e-3
s-1 fitted to these runs, reproduces every run below pH 3.2; I have read that claim and their figures; the engine had
yet to compute any model value for these runs.

**Measured:** the pseudo-first-order constants kobs of De Laat & Gallard 1999, Table 3 (six runs, pH 1.0 to 3.18,
200 uM Fe(III), 10 mM H2O2) and Table 4 (sixteen runs at pH 3.0, Fe(III) 50 to 1000 uM, H2O2 0.2 to 955 mM); 25 °C,
I = 0.1 M (HClO4/NaClO4). Runs 7 to 9 of Table 3 (pH 3.25 and above) stay out: the paper measured iron precipitating
in them (dissolved Fe(III) below the 200 uM added) and limits the model to pH below 3.2.

**Computed:** the mechanism of their Table 2 solved in time with every step written as a mass-action rate, the
hydrolysis and hydroperoxy complexes held at equilibrium, and [H+] = 10^-pH / 0.85 (the activity coefficient the
authors used at I = 0.1 M). kobs is the least-squares slope through the origin of ln([H2O2]/[H2O2]0) against time
over the first 20% of H2O2 removal, the paper's eq 13 and range. The constant k11a (Fe(III) + HO2) is printed as an
upper bound, below 2e3 M-1 s-1; the engine uses 1e3, the value the same group used in 1998 (Gallard et al., New J.
Chem. 22, 263, Table 1), and repeats every comparison at 0 and at 2e3.

**Comparisons:**

1. The published model, k6 = 2.7e-3 s-1, on all 22 runs. These runs set k6, so this comparison confirms the engine
   reproduces the published model.
2. Held out: k6 fitted (least squares in log kobs) to Table 3 alone, then used to predict the 16 runs of Table 4;
   and fitted to Table 4 alone, then used to predict the 6 runs of Table 3. The predicted runs are separate
   experiments, held apart from setting k6.

**Rule** (the project's rule for rates, set in the plan before any data were read): for each compared set,
|log10(predicted/measured)| ≤ 0.3 (a factor of 2) for at least 80% of the runs, and the median |log10| ≤ 0.15.

**Verdict:** passed when comparison 1 and both directions of comparison 2 meet the rule with k11a = 1e3, and again
with k11a = 0 and with k11a = 2e3.

**Also reported:** k6 refitted to all 22 runs, beside the published 2.7e-3.
