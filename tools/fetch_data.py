"""Download a real point-mutant structure series from the RCSB PDB.

Default target is T4 lysozyme (wild-type 2LZM). The Matthews lab solved an
unusually large series of single-point mutants of this protein against a common
wild-type, which makes it one of the few real batches that matches what
StructScan is built to screen.

The script takes the reference entry's sequence from the RCSB data API, runs a
sequence-similarity search to find near-identical entries, and downloads them.

    python tools/fetch_data.py
    python tools/fetch_data.py --reference 1STN --outdir data/staph_nuclease --limit 60

Caveat worth stating in any write-up: these are independent crystal structures.
Differences in crystal packing, resolution, bound ligands and crystallisation
conditions all contribute to the measured deviation alongside the mutation, so
a high-ranking entry is not automatically a mutation effect.

PDB coordinate data is released into the public domain (CC0).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import List, Optional, Sequence

DATA_API = "https://data.rcsb.org/rest/v1/core/polymer_entity/{entry}/{entity}"
SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2/query"
FILE_URL = "https://files.rcsb.org/download/{entry}.pdb"

USER_AGENT = "StructScan/1.0 (educational structural-biology pipeline)"
DEFAULT_REFERENCE = "2LZM"  # T4 lysozyme, wild-type


def _get(url: str, timeout: int = 30) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def _post_json(url: str, payload: dict, timeout: int = 30) -> dict:
    body = json.dumps(payload).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=body,
        headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read())


def fetch_sequence(entry: str, entity: str = "1") -> str:
    """Canonical one-letter sequence of a polymer entity."""
    data = json.loads(_get(DATA_API.format(entry=entry.upper(), entity=entity)))
    sequence = data["entity_poly"]["pdbx_seq_one_letter_code_can"]
    return "".join(sequence.split())


def search_similar(sequence: str, identity: float = 0.9,
                   evalue: float = 0.1, rows: int = 500) -> List[str]:
    """Four-character PDB IDs of entries whose sequence matches at >= identity."""
    payload = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": evalue,
                "identity_cutoff": identity,
                "sequence_type": "protein",
                "value": sequence,
            },
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": rows},
            "results_verbosity": "compact",
            "scoring_strategy": "sequence",
        },
    }
    result = _post_json(SEARCH_API, payload)

    identifiers: List[str] = []
    for item in result.get("result_set", []):
        # "compact" gives plain strings such as "2LZM_1"; verbose gives dicts.
        raw = item if isinstance(item, str) else item.get("identifier", "")
        entry_id = raw.split("_")[0].upper()
        if entry_id and entry_id not in identifiers:
            identifiers.append(entry_id)
    return identifiers


def download_entry(entry: str, destination: Path, timeout: int = 30) -> bool:
    """Fetch one PDB-format file. Returns False if the entry has no .pdb."""
    try:
        payload = _get(FILE_URL.format(entry=entry.upper()), timeout=timeout)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return False  # too large for PDB format, or obsoleted
        raise
    destination.write_bytes(payload)
    return True


def fetch(outdir: Path, reference: str = DEFAULT_REFERENCE, entity: str = "1",
          limit: int = 60, identity: float = 0.9, delay: float = 0.2) -> int:
    outdir = Path(outdir)
    variants_dir = outdir / "variants"
    variants_dir.mkdir(parents=True, exist_ok=True)

    reference = reference.upper()
    print(f"Reference entry : {reference}")

    sequence = fetch_sequence(reference, entity)
    print(f"Sequence        : {len(sequence)} residues")

    reference_path = outdir / "wildtype_reference.pdb"
    if not download_entry(reference, reference_path):
        print(f"error: {reference} is not available in PDB format", file=sys.stderr)
        return 1
    print(f"Wild-type       -> {reference_path}")

    candidates = [e for e in search_similar(sequence, identity=identity) if e != reference]
    print(f"Similar entries : {len(candidates)} at >={identity:.0%} identity"
          f" (downloading up to {limit})")

    downloaded = skipped = 0
    for entry in candidates:
        if downloaded >= limit:
            break
        destination = variants_dir / f"{entry}.pdb"
        if destination.exists():
            downloaded += 1
            continue
        try:
            ok = download_entry(entry, destination)
        except (urllib.error.URLError, TimeoutError) as exc:
            print(f"  {entry}: network error ({exc}); skipped", file=sys.stderr)
            skipped += 1
            continue
        if ok:
            downloaded += 1
            print(f"  [{downloaded:>3}/{limit}] {entry}")
        else:
            skipped += 1
        time.sleep(delay)  # be polite to the RCSB servers

    print(f"\nDownloaded {downloaded} variant structure(s) to {variants_dir}")
    if skipped:
        print(f"Skipped {skipped} entry/entries with no PDB-format file.")
    print(
        "\nRun the analysis with:\n"
        f"  python main.py --reference {reference_path} "
        f"--variants {variants_dir} --trim-outliers"
    )
    return 0 if downloaded else 1


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Download a wild-type reference plus its mutant series from the RCSB PDB.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--reference", default=DEFAULT_REFERENCE,
                        help="4-character PDB ID of the wild-type reference")
    parser.add_argument("--entity", default="1", help="polymer entity ID of the chain of interest")
    parser.add_argument("--outdir", type=Path,
                        default=Path(__file__).resolve().parents[1] / "data" / "t4_lysozyme")
    parser.add_argument("--limit", type=int, default=60, help="maximum variants to download")
    parser.add_argument("--identity", type=float, default=0.9,
                        help="sequence identity cutoff for the similarity search")
    args = parser.parse_args(argv)

    try:
        return fetch(args.outdir, reference=args.reference, entity=args.entity,
                     limit=args.limit, identity=args.identity)
    except urllib.error.URLError as exc:
        print(f"error: could not reach the RCSB PDB ({exc}).", file=sys.stderr)
        print("hint: 'python tools/make_synthetic.py' works offline.", file=sys.stderr)
        return 2
    except KeyError as exc:
        print(f"error: unexpected API response, missing {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
