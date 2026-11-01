# Peroxide engine

Computes what an experimentalist measures when hydrogen peroxide decomposes over a catalyst: the H2O2 left, the rate,
the oxygen formed and dissolved, the gas collected in a burette or the pressure in a sealed vessel, the heat and heat
flow a calorimeter records, and titration volumes. It also applies the fits a lab applies (rate constant, half-life,
initial rate, reaction order, activation energy, Michaelis–Menten constants, turnover).

Every constant traces to a page of a published source, stored with its uncertainty
(`src/peroxide_engine/data/constants/`). Every model reaches the app only after it reproduces measurements under a rule
written down before the comparison (`experiments/<check>/PASS_RULE.md`).

## Status

Phases 0 to 2 are built: thermochemistry from the NBS tables (with the heat of dilution at every strength), water
from the IAPWS equations, oxygen solubility from Benson, Krause and Peterson, the batch vessel (gas burette, sealed or
open, oxygen in balance or crossing the surface at kLa), titrations, the lab's fits, and dissolved iron(III) by the De
Laat and Gallard mechanism (pH 1.0–3.2, 25 °C, perchlorate). Phase 3, the measured catalyst library, is under way.
Seven checks pass and one is missed:

| Check | Result |
|---|---|
| Transcription | 72 NBS values match NIST's digitized tables; every printed equation check value reproduced |
| Heat of decomposition | Solvay's published heats for 35%, 50%, 70% and pure H2O2 within 0.1% |
| Oxygen solubility | Battino's fifteen-source evaluation within 0.44% (5–55 °C); USGS air saturation within 0.07% |
| Exact solutions | H2O2 within 2e-11 of the exact solution; oxygen conserved to 1e-14 |
| Iron model as printed | Every printed hydrolysis constant and ten derived values reproduced (within 3.9%) |
| Iron rates, fitted constant held out | All 22 of De Laat and Gallard's runs within a factor of 2 (largest gap 1.42) |
| Measured catalyst library | ZrO2 and HfO2: every stored value matches its printed text; every related rate, temperature and area the thesis prints is reproduced |
| Iron, another laboratory (missed) | Barb et al.'s rates: 8 of 8 within a factor of 2, typical gap 1.18; their iron(II): typical gap 1.52, above the rule's 1.41, so iron(II) stays off the app |

The catalyst enters as dissolved iron(III), as a measured catalyst from the library (the Catalyst library page lists
every class, with the paper each waiting catalyst needs), or as a first-order rate constant you supply.

Under way: catalysts computed from the crystal alone (`src/peroxide_engine/compute/`). The structure comes from the
Materials Project; pymatgen cuts every stoichiometric, nonpolar, symmetric face; Meta's UMA model relaxes each, picks
the face of lowest surface energy and binds the hydroxyl radical on every site. Measured catalysts are the test set:
the eight oxides Lousada measured ran first, under `experiments/oxide_descriptor_check/PASS_RULE.md`. The faces came
out right; OH binding on the dry face missed its ranking rule (ρ = −0.57 against the measured barriers), so computed
activity stays off the app until a better model passes. The worker reads
one job folder and writes one results file, runs on a GPU when present and resumes after an interruption, so a cluster
runs one folder per GPU:

In the app, Compute a catalyst searches the Materials Project by formula or elements, runs the crystals you tick
through a background queue and shows each face, surface energy and binding site; Data collection gathers every
computed number with the model and run that made it, beside the measured values they are tested against.

From the command line, a screen selects crystals from the Materials Project by filter, as many as match, and splits them into job folders
with a SLURM array script (`run_slurm.sh`, one folder per GPU):

```
uv run pe screen --count --max-elements 3                          # how many oxides match (8,357 at ≤ 40 atoms)
uv run pe screen --max-elements 3 --per-job 50                     # data/screens/O-8357/shard-0000 ... and run_slurm.sh
```

The measured check runs the same way:

```
uv run python experiments/oxide_descriptor_check/prepare.py        # job folders, in this project's Python
TORCHDYNAMO_DISABLE=1 python src/peroxide_engine/compute/worker.py <job folder>   # any Python with fairchem
uv run python experiments/oxide_descriptor_check/check.py           # the verdict
```

Then mechanisms for enzymes, iodide and copper, each after it passes a check against measured rates
(`docs/status.yaml`).

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
