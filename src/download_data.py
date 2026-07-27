"""
Download the raw dataset from Zenodo (record 10.5281/zenodo.15866724) into data/raw/.

Stage 1 reads four Seurat->AnnData `.h5ad` matrices; this module fetches them, plus the
four `*_CoDi_KLD.csv` annotation references that Stage 3 cross-checks against. Standard
library only (`urllib`) so it stays importable without the heavy single-cell stack and
adds no new dependency.

Idempotent: a file already on disk at the correct size is skipped. Each download lands in
a temporary `*.part` file and is renamed only on success, so an interrupted transfer is
never mistaken for a complete file.

    python -m src.download_data           # fetch everything missing (~832 MB h5ad + ~2.4 MB csv)
    python -m src.download_data --no-csv  # only the four .h5ad matrices
    python -m src.download_data --force   # re-download even if already present
"""
from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

from . import config as cfg

ZENODO_RECORD = "15866724"
_BASE = f"https://zenodo.org/api/records/{ZENODO_RECORD}/files/{{name}}/content"

# Exact byte sizes from the Zenodo record -- used both to skip already-complete
# files and to detect a truncated/corrupt download. (Verified 2026-07-27.)
FILE_SIZES = {
    "filtered_feature_bc_matrix.h5ad":                  201_585_802,
    "filtered_feature_bc_matrix_Sample2.h5ad":          272_639_044,
    "filtered_feature_bc_matrix_Sample3.h5ad":          166_968_318,
    "filtered_feature_bc_matrix_Sample4.h5ad":          199_802_137,
    "filtered_feature_bc_matrix_CoDi_KLD.csv":              616_738,
    "filtered_feature_bc_matrix_Sample2_CoDi_KLD.csv":      932_462,
    "filtered_feature_bc_matrix_Sample3_CoDi_KLD.csv":      449_929,
    "filtered_feature_bc_matrix_Sample4_CoDi_KLD.csv":      473_063,
}


def human(n: float) -> str:
    """Bytes -> a compact human-readable string (e.g. 201585802 -> '192.2 MB')."""
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} GB"


def _codi_name(h5ad_name: str) -> str:
    """filtered_feature_bc_matrix_Sample2.h5ad -> ..._Sample2_CoDi_KLD.csv"""
    return h5ad_name[: -len(".h5ad")] + "_CoDi_KLD.csv"


def h5ad_files() -> list[str]:
    """The four raw matrix filenames, taken straight from config.SAMPLES."""
    return list(cfg.SAMPLES.values())


def csv_files() -> list[str]:
    """The four CoDi annotation-reference CSVs that pair with the matrices."""
    return [_codi_name(f) for f in cfg.SAMPLES.values()]


def missing_h5ad() -> list[str]:
    """Raw matrix files not yet present on disk."""
    return [n for n in h5ad_files() if not (cfg.RAW_DIR / n).exists()]


def raw_data_present() -> bool:
    """True when all four raw matrices are on disk at their expected size."""
    for name in h5ad_files():
        path = cfg.RAW_DIR / name
        if not path.exists():
            return False
        want = FILE_SIZES.get(name)
        if want is not None and path.stat().st_size != want:
            return False
    return True


def _download_one(name: str, force: bool = False) -> Path:
    dest = cfg.RAW_DIR / name
    expected = FILE_SIZES.get(name)
    if dest.exists() and not force and (expected is None or dest.stat().st_size == expected):
        print(f"  [skip] {name} ({human(dest.stat().st_size)})")
        return dest

    url = _BASE.format(name=name)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "ANOA-ANO-downloader/1.0"})
    print(f"  [get ] {name}")
    with urllib.request.urlopen(req) as resp:  # noqa: S310 (trusted Zenodo host)
        total = int(resp.headers.get("Content-Length") or expected or 0)
        done = 0
        with open(tmp, "wb") as fh:
            while True:
                buf = resp.read(1 << 20)  # 1 MB chunks
                if not buf:
                    break
                fh.write(buf)
                done += len(buf)
                if total:
                    sys.stdout.write(
                        f"\r         {100 * done / total:5.1f}%  ({human(done)}/{human(total)})"
                    )
                    sys.stdout.flush()
    if total:
        sys.stdout.write("\n")

    if expected is not None and tmp.stat().st_size != expected:
        got = tmp.stat().st_size
        tmp.unlink(missing_ok=True)
        raise OSError(
            f"{name}: downloaded {got} bytes but expected {expected}. "
            "Transfer was likely interrupted -- re-run to retry."
        )
    tmp.replace(dest)  # atomic rename over any previous copy
    return dest


def download_all(include_csv: bool = True, force: bool = False) -> list[Path]:
    """Fetch every raw file (matrices, and CoDi CSVs unless include_csv=False)."""
    names = h5ad_files() + (csv_files() if include_csv else [])
    cfg.RAW_DIR.mkdir(parents=True, exist_ok=True)
    pending = sum(
        FILE_SIZES.get(n, 0) for n in names
        if force or not (cfg.RAW_DIR / n).exists()
    )
    print(f"Zenodo record {ZENODO_RECORD}  ->  {cfg.RAW_DIR}")
    if pending:
        print(f"Downloading ~{human(pending)} ({len(names)} file(s) total, existing ones skipped).")
    else:
        print("All files already present.")
    out = [_download_one(n, force=force) for n in names]
    print(f"Done. {len(out)} file(s) in {cfg.RAW_DIR}.")
    return out


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Download the raw Zenodo dataset into data/raw/.")
    p.add_argument("--no-csv", action="store_true", help="Only the four .h5ad matrices (skip CoDi CSVs).")
    p.add_argument("--force", action="store_true", help="Re-download even if a correct-size file exists.")
    args = p.parse_args(argv)
    try:
        download_all(include_csv=not args.no_csv, force=args.force)
    except Exception as exc:  # noqa: BLE001 -- surface any failure as a clean nonzero exit
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
