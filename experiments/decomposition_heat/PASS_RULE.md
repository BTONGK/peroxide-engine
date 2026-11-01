# Pass rule: heat released by H2O2 solutions

Written 2026-10-31. Before writing it I had computed the 35% row by hand while choosing sources (it agreed), and I had
seen that Solvay's 70% per-mole value (23.4 kcal/mol) differs from what its own 70% per-pound value implies
(23.0 kcal/mol). The rule follows from the stated uncertainties alone.

**Computed:** heat released when a solution of each strength decomposes completely to liquid water and O2 gas at
25 °C, from the NBS dilution series (`thermo.decomposition_heat_of_solution`).

**Compared with:** every heat Solvay prints on page 3 of its Safety & Handling data sheet (HH-2323, 08/2019): the
three per-pound-of-solution heats, the per-pound heat of 100% H2O2, and the three per-mole heats.

**Allowed difference** for each value: half a unit of its last printed digit (rounding), plus the NBS tables'
uncertainty of 0.8 kJ per mol of H2O2 (80 units of the last digit of -191.17 and -187.78, Wagman et al. 1982 §4.2),
expressed in the printed unit. Btu are International Table Btu (2.326 kJ/kg per Btu/lb) and calories thermochemical
(4.184 J).

**Verdict:** passed when all four per-pound heats fall within their allowed difference; these are the values safety
calculations use. The per-mole values are reported beside them with their own result.

**Also measured:** the interpolation error of the dilution series, by leaving each interior NBS point out in turn
and predicting it from the rest. It must stay below 0.1 kJ/mol (an eighth of the tables' uncertainty) for the
interpolation error to stay immaterial.
