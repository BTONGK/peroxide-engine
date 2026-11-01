"""Write the UMA job folders for the oxide check: the eight measured oxides in their measured phases, OH on the face
with the lowest surface energy. Run: uv run python experiments/oxide_descriptor_check/prepare.py
then run each folder in runs/ with compute/worker.py (or check.py --run on this machine)."""

from pathlib import Path

import yaml

from peroxide_engine.compute import jobs

HERE = Path(__file__).parent
data = yaml.safe_load((HERE / "data.yaml").read_text())
MODELS = {"uma-s-1p1_omat": ("~/uma/checkpoints/uma-s-1p1.pt", "omat"),
          "uma-s-1p2p1_oc22": ("~/uma/checkpoints/uma-s-1p2p1.pt", "oc22")}


def main():
    entries = []
    for ox in data["oxides"]:
        e = jobs.material_entry(ox["formula"], ox["space_group_number"])
        print(e["label"], e["id"], "faces:", len(e["slabs"]), "set aside:", e["slabs_set_aside"])
        entries.append(e)
    for name, (ckpt, task) in MODELS.items():
        jobs.write_job(HERE / "runs" / name, entries, ckpt, task, ["OH"])


if __name__ == "__main__":
    main()
