"""The queue of computed runs: each queued run, oldest first, through the UMA worker, one at a time, until the queue is
empty. The app starts it (screens.start_runner); `python -m peroxide_engine.compute.runner` starts it by hand."""

from __future__ import annotations

import os

from peroxide_engine.compute import jobs, screens


def main() -> None:
    pidfile = screens.data_dir() / "runner.pid"
    pidfile.parent.mkdir(parents=True, exist_ok=True)
    pidfile.write_text(str(os.getpid()))
    try:
        while True:
            queued = sorted((r for r in screens.runs() if screens.state(r) == "queued"),
                            key=lambda r: r["meta"]["created"])
            if not queued:
                break
            folder = queued[0]["folder"]
            code = jobs.launch(folder).wait()
            after = screens._status(folder).get("state")
            if after not in ("finished", "stopped"):
                screens._write_status(folder, state="failed" if code else "stopped",
                                      message=f"The worker ended with code {code}; worker.log holds its output")
    finally:
        pidfile.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
