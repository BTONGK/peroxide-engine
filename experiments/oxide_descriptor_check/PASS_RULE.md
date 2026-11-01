# Pass rule: computed OH binding against measured oxide activity

Written 2026-10-31, before any UMA energy for these oxides was computed.

## Question

Does an energy computed from the crystal alone rank how fast oxides decompose H2O2? If it does, the engine can
screen any oxide by formula, and the measured oxides become its test set.

## What was seen first

- The measured values below (Lousada 2013 thesis, Tables 3 and 4).
- The thesis's own DFT trend (Figure 49, p. 98): on small hydroxylated clusters, a more strongly bound pair of HO
  radicals goes with a lower activation energy.
- UMA's error for OH on oxides against OC22's DFT, from catalyst-engine's oxide check: mean 0.21 eV in the
  materials mode (omat) of UMA-S 1.1.

## Computation, fixed in advance

1. **Crystal.** For each oxide, the Materials Project entry in the measured phase's space group with the lowest
   energy above the hull: ZrO2 P2₁/c (14), TiO2 anatase I4₁/amd (141), Y2O3 Ia-3 (206), Fe2O3 R-3c (167),
   CeO2 Fm-3m (225), CuO C2/c (15), HfO2 P2₁/c (14), Gd2O3 Ia-3 (206).
2. **Surface.** Among Miller indices up to 1, every stoichiometric, nonpolar, symmetric slab, relaxed by UMA. The
   face with the lowest surface energy, γ = (E_slab − N·E_bulk per atom) / 2A, is the one a powder mostly shows,
   and the one used.
3. **Binding.** That face repeated to at least 7 Å across, its lower half held at the relaxed positions. OH goes on
   every site the site finder returns: a top site over the highest atom of each element within 2.5 Å of the face
   (so a metal ion under an oxygen layer counts), a bridge and up to two hollows. Each relaxes to 0.05 eV/Å after a
   0.03 Å nudge and is set against the clean face. The descriptor is the lowest value among sites where the OH
   stays whole and bound:
   D = 2 × E_bind(OH) + 2 E(H2O) − E(H2) − E(H2O2), the energy of H2O2 → 2 OH* at wide spacing.
4. **Model.** UMA-S 1.1, materials mode (omat), decides the verdicts. UMA-S 1.2.1's oxide mode (oc22) is reported
   beside it and leaves the verdicts as they are.

## Measured values

- Activation energy Ea, kJ/mol (Table 4, p. 43; 8 oxides): ZrO2 33, TiO2 37, Y2O3 44, Fe2O3 51, CeO2 40, CuO 76,
  HfO2 60, Gd2O3 63.
- Rate per surface area k2 at 298 K, m/s (Table 3, p. 42; 6 oxides, CeO2 and CuO exchanged as Figure 16 and the
  table's own k1 values show): ZrO2 2.39e-10, Fe2O3 3.0e-9, CeO2 1.23e-9, CuO 2.80e-8, HfO2 2.78e-9,
  Gd2O3 9.4e-10.

## Rules

| Verdict | Passes when | Opens in the app |
|---|---|---|
| Ea ranking (the check's verdict) | n = 8, Spearman ρ(D, Ea) ≥ 0.6, exact one-sided p ≤ 0.05 | oxides by formula, ranked by computed D, labeled as an Ea ranking |
| Rate ranking | n = 6, ρ(D, log k2) ≤ −0.6, exact one-sided p ≤ 0.05 | the ranking labeled as a rate ranking |
| Ea value | leave-one-out line Ea = a·D + b: at least 6 of 8 oxides within ±6 kJ/mol | a computed Ea, with its leave-one-out error as the band |

A positive ρ in the first rule means stronger OH binding goes with a lower barrier, the direction Figure 49 shows.
A failed rule stays in the Methodology limits with its numbers, and the next model (hydroxylated faces, the OH pair,
the barrier itself by a nudged elastic band) faces the same rule.
