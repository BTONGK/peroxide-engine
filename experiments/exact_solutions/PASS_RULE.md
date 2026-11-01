# Pass rule: exact solutions and conservation

Written 2026-10-31, before the script first ran.

These are checks of the arithmetic, so the tolerances sit near the solver's precision.

1. **H2O2 against the exact solution.** With a first-order catalyst, the solver's [H2O2] must match c0·exp(-kt)
   within 1e-7 (relative) at every time, in every gas arrangement with finite kLa and with a rate law passed as a
   function.
2. **Oxygen balance.** Starting O2 plus O2 formed must equal O2 in the liquid, in the gas and released, within 1e-9 of
   the total, at every time and in all six arrangements (burette, closed, open; equilibrium or finite kLa).
3. **Equilibrium as the limit of fast transfer.** With kLa at 1000 times k, gas readouts must agree with the
   equilibrium arrangement within 0.1%.
4. **Heat bookkeeping.** The heat over a run must equal (H2O2 decomposed × 100.67 kJ/mol, the enthalpy to dissolved
   O2) minus (O2 moved to the gas × its enthalpy of desorption, plus water evaporated into the new gas at the burette),
   within 1e-4.
5. **The lab's fits recover known answers.** On noise-free computed data, the first-order fit returns k, the order
   fit returns 1 and 2 for first- and second-order series, the Arrhenius fit returns the Ea used, and the
   Michaelis–Menten fit returns Vmax and Km, each within 1e-6 (relative).

**Verdict:** passed when all five hold.
