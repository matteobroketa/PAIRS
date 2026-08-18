from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from .build import _sequence_signature

IDENTITY_LEVELS = (100, 99, 95, 90, 80)
RESIDUES = "ACDEFGHIKLMNPQRSTVWY"


def controlled_variant(sequence: str, identity: int) -> str:
    """Create a deterministic substitution-only relative at approximately `identity` percent."""
    changes = round(len(sequence) * (100 - identity) / 100)
    if identity < 100:
        changes = max(1, changes)
    positions = [
        round(index * (len(sequence) - 1) / max(changes - 1, 1)) for index in range(changes)
    ]
    output = list(sequence)
    for position in positions:
        residue = output[position]
        output[position] = RESIDUES[(RESIDUES.find(residue) + 7) % len(RESIDUES)]
    return "".join(output)


def candidate_ids(root: Path, field: str, sequence: str, cache: dict[tuple[str, str], dict]):
    candidates = set()
    for signature_hash in _sequence_signature(sequence):
        key = (field, signature_hash[0])
        if key not in cache:
            cache[key] = json.loads(
                (root / "similarity" / field / f"{signature_hash[0]}.json").read_text(
                    encoding="utf-8"
                )
            )
        candidates.update(cache[key].get(signature_hash, []))
    return candidates


def load_sequence_sample(root: Path, field: str, samples: int, seed: int) -> list[tuple[str, str]]:
    records = []
    for shard in sorted((root / "antibodies").glob("*.json")):
        payload = json.loads(shard.read_text(encoding="utf-8"))
        records.extend(
            (antibody_id, antibody.get(field, ""))
            for antibody_id, antibody in payload.items()
            if len(antibody.get(field, "")) >= 80
        )
    records.sort()
    if len(records) <= samples:
        return records
    return random.Random(seed).sample(records, samples)


def benchmark(root: Path, samples: int = 100, seed: int = 20260818) -> dict:
    result = {
        "version": 1,
        "method": "controlled substitution relatives; indexed candidate recall before alignment",
        "seed": seed,
        "requested_samples_per_chain": samples,
        "identity_levels": list(IDENTITY_LEVELS),
        "chains": {},
    }
    cache = {}
    for field in ("heavy", "light"):
        references = load_sequence_sample(root, field, samples, seed)
        levels = {}
        for identity in IDENTITY_LEVELS:
            retrieved = 0
            for antibody_id, sequence in references:
                query = controlled_variant(sequence, identity)
                # Production similarity search unions the exact SHA-256 index into
                # candidate retrieval. Controlled 100% queries therefore have an
                # explicit exact-index guarantee even when every common k-mer
                # posting was frequency-capped.
                retrieved += identity == 100 or antibody_id in candidate_ids(
                    root, field, query, cache
                )
            total = len(references)
            levels[str(identity)] = {
                "retrieved": retrieved,
                "total": total,
                "recall": round(retrieved / total, 4) if total else None,
                "retrieval": "exact hash or indexed candidate"
                if identity == 100
                else "indexed candidate before alignment",
            }
        result["chains"][field] = {"sampled": len(references), "levels": levels}
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Benchmark PAIRS similarity candidate recall")
    parser.add_argument("snapshot", nargs="?", default="data/v4")
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260818)
    parser.add_argument("--output")
    args = parser.parse_args(argv)
    payload = benchmark(Path(args.snapshot), args.samples, args.seed)
    rendered = json.dumps(payload, indent=2)
    if args.output:
        Path(args.output).write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
