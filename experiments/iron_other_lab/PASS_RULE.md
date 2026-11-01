# Pass rule: iron(III) rates and ferrous iron measured by another laboratory

Written 2026-10-31, before the iron code existed. While reading Barb et al. I estimated two of their conditions by
hand with the simplified steady-state formula of De Laat & Gallard (eqs 20 to 22, with the factor 2 of eq 21): kd at
[H+] = 0.0049 M came out about 25% above Barb's values and ferrous iron at 4.85 mM Fe(III) about 30% below theirs.
The rule below is the project's rule for rates, set in the plan before these estimates.

**Measured** (Barb, Baxendale, George & Hargrave 1951, Trans. Faraday Soc. 47, 591; 25 °C):

- Table I: the second-order constant kd for four solutions ([H+] 0.204 to 0.0049 M, Fe(III) 19.6 to 0.98 mM),
  four values from oxygen evolution at 0.75 M H2O2 and four from permanganate titration at 1.4 M H2O2.
- Table IV: the steady ferrous iron [Fe2+]s in 16 mixtures at [H+] = 0.025 M, Fe(III) 4.85 to 14.55 mM and H2O2
  0.025 to 0.15 M.

Barb's solutions differ from the model's calibration: nitrate in place of perchlorate, ionic strength 0.435 M in
place of 0.1 M, peroxide up to 1.4 M and iron up to 19.6 mM. The engine uses the model exactly as published for 0.1 M
perchlorate, with every constant unchanged and k6 = 2.7e-3 s-1. This comparison shows how far the model carries
beyond the solutions it was built on.

**Computed:**

- [H+]: Barb lists the acid added. The engine solves the free [H+] at the start from the proton balance (each
  hydroxo and hydroperoxo iron species releases its protons) and holds it fixed through the run, as the model holds
  pH fixed.
- Oxygen evolution (Table Ia): Barb measured rates at nearly constant H2O2. The engine reports
  kd = -(d[H2O2]/dt) / ([H2O2][Fe(III)]0) at the moment ferrous iron reaches its highest value.
- Titration (Table Ib): kd = kobs / [Fe(III)]0, with kobs the least-squares slope of ln[H2O2] against time over the
  first half of the reaction (Barb plotted first-order lines and leaves their range unstated). The slope over the
  first 20% is reported beside it.
- Ferrous iron (Table IV): [Fe2+] three minutes after mixing, starting from 1e-5 M Fe(II) as Barb did ("Generally
  about 10^-5 M ferrous ion was added", p. 608; the mixture "was allowed to stand for a few minutes"). Values at one
  and ten minutes are reported beside it.
- k11a = 1e3 M-1 s-1, repeated at 0 and 2e3, as in the same-laboratory check.

**Rule:** for the 8 kd values, and separately for the 16 ferrous iron values, |log10(predicted/measured)| ≤ 0.3 for
at least 80% of the values and the median |log10| ≤ 0.15.

**Verdict:** passed when both sets meet the rule with k11a = 1e3, and again with k11a = 0 and with k11a = 2e3.

**Use in the app:** the app offers the iron model only inside the conditions of De Laat & Gallard's runs, where the
same-laboratory check tests it. This check is reported on the Model checks page with its verdict, as evidence of
how the model fares in another laboratory's solutions.
