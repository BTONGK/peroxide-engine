# Peroxide engine

Computes what an experimentalist measures when hydrogen peroxide decomposes over a catalyst: the H2O2 left, the rate,
the oxygen formed and dissolved, the gas collected in a burette or the pressure in a sealed vessel, the heat and heat
flow a calorimeter records, and titration volumes. It also applies the fits a lab applies (rate constant, half-life,
initial rate, reaction order, activation energy, Michaelis–Menten constants, turnover).

Every constant traces to a page of a published source, stored with its uncertainty
(`src/peroxide_engine/data/constants/`). Every model reaches the app only after it reproduces measurements under a rule
written down before the comparison (`experiments/<check>/PASS_RULE.md`).

## Status

Phases 0 and 1 are built: thermochemistry from the NBS tables (with the heat of dilution at every strength), water
from the IAPWS equations, oxygen solubility from Benson, Krause and Peterson, the batch vessel (gas burette, sealed or
open, oxygen in balance or crossing the surface at kLa), titrations and the lab's fits. Four checks pass:

| Check | Result |
|---|---|
| Transcription | 72 NBS values match NIST's digitized tables; every printed equation check value reproduced |
| Heat of decomposition | Solvay's published heats for 35%, 50%, 70% and pure H2O2 within 0.1% |
| Oxygen solubility | Battino's fifteen-source evaluation within 0.44% (5–55 °C); USGS air saturation within 0.07% |
| Exact solutions | H2O2 within 2e-11 of the exact solution; oxygen conserved to 1e-14 |

Today the catalyst enters as a first-order rate constant you supply. Mechanisms for dissolved iron, enzymes, iodide,
copper, oxides and metals follow, each after it passes a check against measured rates (`docs/status.yaml`).

## Use

```bash
uv sync
uv run pe ui            # the app, on the first free port from 8502 (or --port)
uv run pe checks        # rerun every model check
uv run pytest
```

```python
from peroxide_engine import vessel
run = vessel.simulate(vessel.Setup(c0_M=0.1, k_per_s=1e-3, headspace="burette"))
run.readouts            # each number with its units and how it was computed
run.series              # the time series
```

## Papers

Published papers stay in the git-ignored `cache/`; the repository holds the digitized numbers, the scripts and the
citations.
