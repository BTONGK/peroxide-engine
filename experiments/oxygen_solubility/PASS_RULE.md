# Pass rule: oxygen solubility

Written 2026-10-31, after reading Battino et al.'s Table 1 and before running this script.

**Computed:** the mole fraction of O2 dissolved in water at 101.325 kPa O2 from Benson, Krause and Peterson (1979)
(`oxygen.mole_fraction_1atm`), and from it the O2 in water saturated with air (`oxygen.air_saturation`).

**Comparison 1, independent evaluation:** Battino, Rettich and Tominaga's recommended values (1983, Table 1,
column b), smoothed from 149 measurements in fifteen sources. Allowed difference: 0.5% at every tabulated
temperature from 5 to 55 °C. That range is the app's temperature range. I set it after reading Table 1: above 55 °C the
two evaluations part by 0.7–2.7%, so the app stops at 55 °C.

**Comparison 2, air saturation:** Benson and Krause's 1984 formula for O2 in fresh water under moist air at 1 atm
(mg/L), as adopted by the USGS (Technical Memorandum 2011.03, eq. 7). Allowed difference: 0.1% from 0.01 °C (the IAPWS water equations start at the triple point) to
40 °C, the formula's range, in steps of 5 °C. Both come from the same laboratory's measurements, so this checks the engine's conversion from
mole fraction to mg/L (air's O2 fraction, water vapor, density, non-ideality).

**Verdict:** passed when both comparisons hold at every temperature.

**Reported beside it:** the enthalpy of solution at 25 °C from the temperature slope, against Battino et al.'s
Table 2 (-12.11 kJ/mol) and the NBS tables' O2(ao) (-11.7 kJ/mol, uncertain by up to 8 kJ/mol under the tables'
rule).
