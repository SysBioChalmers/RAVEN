#!/usr/bin/env python
"""Populate RAVEN's software/ with every platform's external binaries.

Downloads the per-platform binary ZIPs from toolbox-assets and lays them out under
``<raven>/software/<tool>/`` using RAVEN's suffix scheme — bare name = Linux,
``.mac`` = macOS, ``.exe`` = Windows. With all platforms present in one tree, the
whole RAVEN checkout can be zipped into a single **OS-independent** offline
distribution (see the release-bundle workflow).

Run from the RAVEN root, or pass ``--raven <dir>``. Requires network access to
https://github.com/CompMetEng/toolbox-assets .
"""
from __future__ import annotations

import argparse
import io
import urllib.request
import zipfile
from pathlib import Path

TOOLBOX_ASSETS = "https://github.com/CompMetEng/toolbox-assets/releases/download"
# toolbox-assets platform key -> RAVEN binary suffix.
PLATFORMS = {"linux-x86_64": "", "macos-arm64": ".mac", "windows-x86_64": ".exe"}

# (software subdir, toolbox-assets bundle, {platform|"_": version}, {bare executable names})
CLI_TOOLS = [
    ("blast+",  "blast",   {"_": "2.17.0"},                          {"blastp", "makeblastdb"}),
    ("diamond", "diamond", {"_": "2.1.17"},                          {"diamond"}),
    ("hmmer",   "hmmer",   {"windows-x86_64": "3.3.2", "_": "3.4.0"}, {"hmmsearch"}),
]


def _fetch(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as resp:  # noqa: S310 (pinned toolbox-assets URLs)
        return resp.read()


def _extract(zip_bytes: bytes, dest: Path, rename: dict[str, str]) -> None:
    """Write each file member of the ZIP into ``dest``, applying ``rename`` and modes."""
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
        for info in zf.infolist():
            if info.is_dir():
                continue
            target = dest / rename.get(info.filename, info.filename)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(zf.read(info))
            mode = (info.external_attr >> 16) & 0o777
            if mode:
                target.chmod(mode)  # no-op for execute bits on Windows; honoured on the Linux runner


def assemble(raven: Path) -> None:
    soft = raven / "software"
    for subdir, bundle, vmap, exes in CLI_TOOLS:
        dest = soft / subdir
        for plat, suffix in PLATFORMS.items():
            version = vmap.get(plat, vmap["_"])
            asset = f"{bundle}-{version}-{plat}.zip"
            print(f"  {subdir}: {asset}")
            # macOS ships bare names; rename them to <name>.mac. Linux keeps bare names,
            # Windows already carries .exe — so only the macOS members need renaming.
            rename = {e: e + ".mac" for e in exes} if suffix == ".mac" else {}
            _extract(_fetch(f"{TOOLBOX_ASSETS}/{bundle}-{version}/{asset}"), dest, rename)


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--raven", type=Path, default=Path.cwd(), help="RAVEN root (default: cwd)")
    args = p.parse_args(argv)
    assemble(args.raven)
    print("Done.")


if __name__ == "__main__":
    main()
