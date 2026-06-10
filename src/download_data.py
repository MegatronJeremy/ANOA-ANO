"""
Download the dataset from Zenodo into data/raw/.

Usage:  python src/download_data.py

Zenodo record: https://zenodo.org/records/15866724  (DOI 10.5281/zenodo.15866724)

NOTE: this fetches the record's file manifest from the Zenodo API and downloads
every file. After it runs, open data/raw/, see what the .h5ad files are actually
called, and update the SAMPLES dict in src/config.py to match.
"""
import urllib.request
import json
from pathlib import Path

RECORD_ID = "15866724"
API_URL = f"https://zenodo.org/api/records/{RECORD_ID}"
RAW_DIR = Path(__file__).resolve().parents[1] / "data" / "raw"


def main():
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(API_URL) as r:
        record = json.load(r)

    files = record.get("files", [])
    print(f"{len(files)} files in record {RECORD_ID}:")
    for f in files:
        name = f["key"]
        url = f["links"]["self"]
        dest = RAW_DIR / name
        if dest.exists():
            print(f"  [skip] {name} already present")
            continue
        print(f"  [get ] {name} ...")
        urllib.request.urlretrieve(url, dest)
    print(f"\nDone. Files in {RAW_DIR}. Now update SAMPLES in src/config.py.")


if __name__ == "__main__":
    main()
