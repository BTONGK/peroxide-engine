# Pass rule: each measured catalyst in the library holds together with its own source

Written 2026-10-31, before the library code and this check's script existed. While reading Lousada's thesis to
decide which oxides it supports, I worked the ZrO2 and HfO2 arithmetic below by hand; that reading is also how the
CeO2, CuO, Fe2O3, Gd2O3, TiO2 and Y2O3 entries came to be held for the journal paper. The rule below is fixed now.

**What the library holds:** for each catalyst marked ready, the rate law a laboratory measured, the values as printed,
and the conditions they were measured over (H2O2, temperature, catalyst amount). The app runs a catalyst only inside
those conditions.

**Checks, for every ready catalyst:**

1. Transcription: each stored value equals the number its printed text states (value and printed uncertainty).
2. The measured rate at a stated amount: where the source prints k1 for a stated surface area, k2 × Sa/V + b computed
   from the stored values must lie within the sum of three allowances: k2's printed uncertainty times Sa/V, half a
   unit in the last printed digit of b, and k1's printed uncertainty.
3. Temperature: where the source prints k1 at two temperatures from one series, the stored activation energy must
   give their ratio within the effect of its printed uncertainty, plus the printed uncertainties of both k1 values.
4. Pre-exponential factor: A × exp(−Ea/RT) at 298.15 K must match k1 at the source's standard amount within the
   effect of Ea's printed uncertainty and half a unit in A's last printed digit.
5. Surface area: the stored specific surface area times a stated mass gives the stated surface area within the
   specific surface area's printed uncertainty.

**Verdict:** passed when every ready catalyst meets every check that its source supports.
