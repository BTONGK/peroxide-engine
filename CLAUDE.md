# peroxide-engine

A computational tool for hydrogen peroxide decomposition only: it computes what an experimentalist measures (rate,
oxygen, pressure, heat, titrations, the lab's fits) for any catalyst class. The user ranks accuracy above features.

- **Every number has a source.** Constants, tables and equations live in `src/peroxide_engine/data/constants/*.yaml`
  with a DOI or URL, a locator (table, page, equation) and an uncertainty with how it was obtained; `sources.py`
  refuses a record without them. Store table values as printed strings (the last digit carries the uncertainty).
  Read each value from the source itself (page images for scans), and record where in `transcribed_from` or the
  locator.
- **Checked models only.** A model joins the app only after it reproduces measurements it was built apart from, under
  a rule in `experiments/<check>/PASS_RULE.md` written before the comparison (say so honestly when something was seen
  first). Each check folder holds `check.py`, `results/summary.json` (written by the script, with `verdict`) and
  `report.yaml` (plain words; `{key}` fills from the summary). `tests/test_checks.py` holds each verdict. What has yet
  to pass goes in `docs/methodology.yaml` under `limits`.
- **Copy is affirmative**: no negation words in app strings, docs yaml, reports, PASS_RULE files or the README;
  `tests/test_copy.py` enforces it.
- **Papers stay local**: PDFs and downloaded tables go in the git-ignored `cache/` (`cache/papers/`, `cache/nbs/`);
  `tests/test_sources.py` fails if anything under `cache/` or any PDF is tracked. The repository is public.
- Update `docs/status.yaml` (the Project status page) at the end of every piece of work.
- The app: `uv run pe ui` (first free port from 8502, so it runs beside catalyst-engine on 8501); pages in `src/peroxide_engine/ui/views/`, registered in `ui/pages.py`. Use a private port
  (e.g. 8778) for checks and stop the server afterwards.
- Catalyst modules (iron, enzymes, iodide, copper, oxides, metals) plug into `vessel.simulate(setup, rate=...)` as rate
  laws or mechanisms; solids will reuse `/Users/bt/catalyst-engine` (Materials Project lookup, UMA worker, ranking
  statistics). The Materials Project key lives in `~/.config/.pmgrc.yaml`; keep it out of the repository and logs.
- Commit only when asked. Use they/them for people.
