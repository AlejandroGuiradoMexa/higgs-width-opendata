"""Download the exactly4lep skim files listed in hwidth.samples to data/raw/."""
import pathlib
import subprocess
import sys

import atlasopenmagic as atom

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from hwidth.samples import RELEASE, SAMPLES, SKIM  # noqa: E402

RAW = pathlib.Path(__file__).resolve().parents[1] / "data" / "raw"


def main():
    atom.set_release(RELEASE)
    RAW.mkdir(parents=True, exist_ok=True)
    for group, keys in SAMPLES.items():
        for key in keys:
            for url in atom.get_urls(str(key), skim=SKIM, protocol="https"):
                url = url.removeprefix("simplecache::")
                out = RAW / url.rsplit("/", 1)[1]
                if out.exists():
                    continue
                print(f"[{group}] {out.name}")
                subprocess.run(["curl", "-sfL", "-o", str(out), url], check=True)


if __name__ == "__main__":
    main()
