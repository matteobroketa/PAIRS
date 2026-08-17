from __future__ import annotations

import json
import sys
from pathlib import Path

EXPECTED_SCHEMA = 3


def _read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def validate(data_dir: Path) -> list[str]:
    errors: list[str] = []
    for name in ["manifest.json", "targets.json"]:
        if not (data_dir / name).exists():
            errors.append(f"missing {name}")
    if errors:
        return errors

    manifest = _read_json(data_dir / "manifest.json")
    targets = _read_json(data_dir / "targets.json")
    antibody_shards = list((data_dir / "antibodies").glob("*.json"))

    if manifest.get("schema_version") != EXPECTED_SCHEMA:
        errors.append("unexpected schema_version")
    if manifest.get("stats", {}).get("targets") != len(targets):
        errors.append("target count mismatch")
    if not antibody_shards:
        errors.append("no antibody shards")

    antibody_ids: set[str] = set()
    antibodies_by_id: dict[str, dict] = {}
    for shard_path in antibody_shards:
        payload = _read_json(shard_path)
        if not isinstance(payload, dict):
            errors.append(f"antibody shard is not an object: {shard_path.name}")
            continue
        for antibody_id in payload:
            if antibody_id in antibody_ids:
                errors.append(f"duplicate antibody id {antibody_id}")
            antibody_ids.add(antibody_id)
            antibodies_by_id[antibody_id] = payload[antibody_id]
            expected_shard = antibody_id[3:5] if antibody_id.startswith("ab_") else None
            if expected_shard and shard_path.stem != expected_shard:
                errors.append(
                    f"antibody {antibody_id} is in shard {shard_path.stem}, expected {expected_shard}"
                )

    if manifest.get("stats", {}).get("antibodies") != len(antibody_ids):
        errors.append("antibody count mismatch")

    target_ids: set[str] = set()
    referenced_antibodies: set[str] = set()
    counted_results = 0
    for target in targets:
        target_id = target.get("id", "")
        if target_id in target_ids:
            errors.append(f"duplicate target id {target_id}")
        target_ids.add(target_id)
        target_dir = data_dir / "targets" / target.get("dir", "")
        index_path = target_dir / "index.json"
        if not index_path.exists():
            errors.append(f"missing target index for {target_id}")
            continue
        target_manifest = _read_json(index_path)
        pages = target_manifest.get("pages", [])
        if len(pages) != target.get("page_count"):
            errors.append(f"target page count mismatch for {target_id}")
        target_result_count = 0
        for page in pages:
            page_path = target_dir / page.get("file", "")
            if not page_path.exists():
                errors.append(f"missing target page {target_id}/{page.get('file')}")
                continue
            rows = _read_json(page_path)
            if len(rows) != page.get("count"):
                errors.append(f"page row count mismatch for {target_id}/{page.get('file')}")
            target_result_count += len(rows)
            for row in rows:
                antibody_id = row.get("antibody", {}).get("id", "")
                if antibody_id not in antibody_ids:
                    errors.append(
                        f"target {target_id} references missing antibody {antibody_id or '<empty>'}"
                    )
                else:
                    referenced_antibodies.add(antibody_id)
                for interaction in row.get("interactions", []):
                    if interaction.get("target_id") != target_id:
                        errors.append(
                            f"interaction {interaction.get('id')} has wrong target_id in {target_id}"
                        )
        if target_result_count != target.get("result_count"):
            errors.append(f"target result count mismatch for {target_id}")
        counted_results += target_result_count

    if not any((data_dir / "antibody-search").glob("*.json")):
        errors.append("no antibody search shards")
    if not any((data_dir / "sequence-search").glob("*.json")):
        errors.append("no sequence search shards")

    for sequence_path in (data_dir / "sequence-search").glob("*.json"):
        payload = _read_json(sequence_path)
        for sequence_hash, matches in payload.items():
            if not sequence_hash.startswith(sequence_path.stem):
                errors.append(f"sequence hash in wrong shard: {sequence_hash}")
            for match in matches:
                antibody_id = match.get("id", "")
                if antibody_id not in antibody_ids:
                    errors.append(
                        f"sequence index {sequence_path.name} references missing antibody {antibody_id}"
                    )

    for field in ("cdrh3", "cdrl3"):
        directory = data_dir / "sequence" / field
        if not directory.exists():
            errors.append(f"missing {field} index directory")
            continue
        for bucket_path in directory.glob("*.json"):
            try: bucket_length = int(bucket_path.stem)
            except ValueError:
                errors.append(f"invalid {field} bucket name: {bucket_path.name}"); continue
            for record in _read_json(bucket_path):
                sequence = record.get("sequence", "")
                if len(sequence) != bucket_length: errors.append(f"{field} sequence in wrong bucket: {sequence}")
                for antibody_id in record.get("antibody_ids", []):
                    if antibody_id not in antibody_ids: errors.append(f"{field} index references missing antibody {antibody_id}")
                    elif sequence != antibodies_by_id[antibody_id].get(field, ""): errors.append(f"{field} index sequence not present on antibody {antibody_id}")

    stats = manifest.get("stats", {})
    if stats.get("interactions", 0) <= 0:
        errors.append("interaction count is zero")
    if targets and not referenced_antibodies:
        errors.append("target pages reference no antibodies")
    if counted_results <= 0:
        errors.append("target result pages are empty")
    return errors


def main() -> int:
    data_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "data/v3")
    errors = validate(data_dir)
    if errors:
        print("\n".join("ERROR: " + error for error in errors), file=sys.stderr)
        return 1
    print(f"Validation OK: {data_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
