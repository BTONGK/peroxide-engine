"""App copy states what is true, with no negation words: every string in the yaml the app shows, the README, and
string literals (docstrings aside) in the package."""

import ast
import re
from pathlib import Path

import yaml

from peroxide_engine import checks

WORDS = re.compile(r"\b(not|no|never|without|none|cannot|nor|neither|nothing|nobody|nowhere|\w+n't)\b", re.I)
ROOT = checks.ROOT


def yaml_strings(node, where=""):
    if isinstance(node, dict):
        for k, v in node.items():
            yield from yaml_strings(v, f"{where}.{k}")
    elif isinstance(node, list):
        for i, v in enumerate(node):
            yield from yaml_strings(v, f"{where}[{i}]")
    elif isinstance(node, str):
        yield where, node


def py_strings(path: Path):
    tree = ast.parse(path.read_text())
    docs = set()
    for n in ast.walk(tree):
        if isinstance(n, (ast.Module, ast.FunctionDef, ast.ClassDef, ast.AsyncFunctionDef)) and n.body:
            first = n.body[0]
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant):
                docs.add(id(first.value))
    for n in ast.walk(tree):
        if isinstance(n, ast.Constant) and isinstance(n.value, str) and id(n) not in docs:
            yield f"line {n.lineno}", n.value


def hits():
    files = [*(ROOT / "docs").glob("*.yaml"), *(ROOT / "experiments").glob("*/report.yaml"),
             *(ROOT / "src" / "peroxide_engine" / "data" / "constants").glob("*.yaml"),
             *(ROOT / "src" / "peroxide_engine" / "data" / "library").glob("*.yaml")]
    out = []
    for f in files:
        for where, text in yaml_strings(yaml.safe_load(f.read_text())):
            out += [f"{f.name}{where}: {m.group(0)}" for m in WORDS.finditer(text) if not m.group(0).isupper()]
    for f in (ROOT / "src" / "peroxide_engine").rglob("*.py"):
        for where, text in py_strings(f):
            out += [f"{f.name} {where}: {m.group(0)}" for m in WORDS.finditer(text) if not m.group(0).isupper()]
    for f in [ROOT / "README.md", *(ROOT / "experiments").glob("*/PASS_RULE.md")]:
        out += [f"{f.name}: {m.group(0)}" for m in WORDS.finditer(f.read_text()) if not m.group(0).isupper()]
    return out


def test_copy_is_affirmative():
    assert hits() == []
