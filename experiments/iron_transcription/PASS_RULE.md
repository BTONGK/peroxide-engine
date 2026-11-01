# Pass rule: the iron model's constants and equations, as printed

Written 2026-10-31, before the iron code existed. While reading the papers I worked several of these values by hand
to learn the papers' conventions: the conversion of pH to [H+] (activity coefficient 0.85 at I = 0.1 M), the
complexed fraction at pH 3, and the rate constants k10app and k11app at pH 3 all matched by hand; the acid constant
linking the two hydroperoxy complexes came out at 1.87e-4 against the printed 1.8e-4. The rule below covers every
value, including that one.

**Computed by the engine** from the constants it stores (De Laat & Gallard 1999, Table 2; Gallard, De Laat & Legube
1999, Table 1):

- Group A, the ionic-strength equations of Gallard et al. Table 1: K1,1 and K1,2 at I = 0.1 M (printed in the same
  table: 2.90e-3 and 7.62e-7) and at I = 0.2 M (printed in the caption of De Laat & Gallard's Figure 4: 2.44e-3 and
  6.47e-7).
- Group B, values the papers derive and print: the complex formation constants written with HO2- (Gallard et al.
  eqs 23 and 24: 1.74e9 and 1.12e8) and the acid constant linking the two complexes (eq 25: 1.8e-4); the share of
  Fe(III) held as hydroperoxy complexes at pH 3.0 with 0.5 M H2O2 and 0.5 mM Fe(III) (27%, Gallard et al. p. 2934);
  the initiation constant ki of De Laat & Gallard eq 20 at pH 1.0, 2.0 and 3.0 with 0.1 M H2O2 (6.9e-5, 5.6e-4,
  1.9e-3) and at pH 3.0 with 1 M H2O2 (1.2e-3); and k10app and k11app at pH 3 (1.3e6 and 6.6e5, p. 2731).

**Allowed difference:** Group A, half a unit of the last printed digit (these are the model's inputs and the
equations must give them exactly). Group B, 5% of the printed value: each is printed to two or three significant
figures after intermediate rounding by the authors, which can move a two-figure value by up to 5%.

**Amended 2026-10-31, before the check script existed:** after writing this rule I worked Group A by hand from the
page image of Table 1. K1,2 at I = 0.1 M comes to 7.6254e-7, which rounds to 7.63; the printed 7.62 equals
10^-6.118, the value of log K1,2 rounded to three decimals, and 2.90e-3 likewise equals 10^-2.538. The authors
rounded log K first. Group A's allowed difference is therefore half a unit in the third decimal of log K (0.0005,
or 0.12% in K), the precision the printed values carry.

**Amended a second time, 2026-10-31, after the first run of the script:** that first run reported K1,1 at
I = 0.2 M as missed, by 0.00007 in log K beyond the 0.0005 allowed. The tolerance above allowed for the authors'
rounding of log K and overlooked the second rounding, of K itself to three significant figures (2.4434e-3 printed
as 2.44e-3). Group A now repeats the authors' two roundings exactly: log K rounded to three decimals, then K rounded
to the printed number of significant figures, must equal the printed value. A coefficient misread by one digit moves
log K by 0.001 or more and still fails this test. The summary keeps the first run's result.

**Verdict:** passed when every Group A and Group B value lies within its allowed difference.

**Also reported:** the second-order constant kd at pH 3 from De Laat & Gallard's eqs 21 and 22 combined, and from
their eq 24 as printed, beside their measured 0.47 ± 0.03 M-1 s-1.
