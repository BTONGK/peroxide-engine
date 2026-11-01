"""Model checks: each folder in experiments/ holds PASS_RULE.md, a script, results/summary.json and report.yaml.

The script writes the summary; the app renders the report with the summary's numbers. A check's verdict is "passed"
or "missed", set by the script from the rule in PASS_RULE.md.
"""

from __future__ import annotations

import json
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
EXPERIMENTS = ROOT / "experiments"


def write_summary(folder: Path, summary: dict) -> Path:
    out = folder / "results" / "summary.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2, default=float) + "\n")
    return out


def load(folder: Path) -> dict:
    report = yaml.safe_load((folder / "report.yaml").read_text())
    summary = json.loads((folder / "results" / "summary.json").read_text())
    return {"folder": folder.name, "report": report, "summary": summary}


def all_checks() -> list[dict]:
    out = [load(f) for f in sorted(EXPERIMENTS.iterdir()) if (f / "report.yaml").exists()]
    return sorted(out, key=lambda c: c["report"].get("order", 99))
