from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
from pathlib import Path

from .build import _sequence_signature
from .model import sequence
from .targets import key

SUPPORTED_SCHEMA = 4


class Snapshot:
    def __init__(self, root: Path):
        self.root = root
        self.manifest = self.read("manifest.json")
        if self.manifest.get("schema_version") != SUPPORTED_SCHEMA:
            raise ValueError(
                f"unsupported snapshot schema {self.manifest.get('schema_version')!r}; "
                f"expected {SUPPORTED_SCHEMA}"
            )
        self.targets = self.read("targets.json")
        self.targets_by_id = {target["id"]: target for target in self.targets}

    def read(self, relative: str | Path):
        path = self.root / relative
        if not path.is_file():
            raise FileNotFoundError(f"snapshot file is missing: {path}")
        return json.loads(path.read_text(encoding="utf-8"))

    def antibody(self, antibody_id: str) -> dict | None:
        if not antibody_id.startswith("ab_") or len(antibody_id) < 5:
            return None
        path = self.root / "antibodies" / f"{antibody_id[3:5]}.json"
        if not path.is_file():
            return None
        return json.loads(path.read_text(encoding="utf-8")).get(antibody_id)

    def resolve_target(self, query: str) -> dict | None:
        if query in self.targets_by_id:
            return self.targets_by_id[query]
        normalized = key(query)
        matches = [
            target
            for target in self.targets
            if normalized
            in {key(target.get("name", "")), *(key(a) for a in target.get("aliases", []))}
        ]
        return matches[0] if len(matches) == 1 else None

    def target_rows(self, target: dict, family: str = "positive") -> list[dict]:
        prefix, count_field = {
            "positive": ("page", "page_count"),
            "functional": ("functional-page", "functional_page_count"),
            "negative": ("negative-page", "negative_page_count"),
        }[family]
        rows = []
        for number in range(1, target.get(count_field, 0) + 1):
            rows.extend(self.read(Path("targets") / target["dir"] / f"{prefix}-{number:03d}.json"))
        return rows

    def descendants(self, target: dict) -> list[dict]:
        output, pending, seen = [], list(target.get("children", [])), set()
        while pending:
            target_id = pending.pop(0)
            if target_id in seen:
                continue
            seen.add(target_id)
            child = self.targets_by_id.get(target_id)
            if child:
                output.append(child)
                pending.extend(child.get("children", []))
        return output

    def exact_sequence(self, amino_acids: str, field: str) -> list[dict]:
        digest = hashlib.sha256(amino_acids.encode()).hexdigest()
        path = self.root / "sequence-search" / f"{digest[:2]}.json"
        if not path.is_file():
            return []
        payload = json.loads(path.read_text(encoding="utf-8"))
        return [match for match in payload.get(digest, []) if match.get("field") == field]

    def similarity(self, amino_acids: str, field: str) -> list[dict]:
        candidate_ids: set[str] = set()
        for signature in _sequence_signature(amino_acids):
            payload = self.read(Path("similarity") / field / f"{signature[0]}.json")
            candidate_ids.update(payload.get(signature, []))
        results = []
        for antibody_id in candidate_ids:
            antibody = self.antibody(antibody_id)
            candidate = antibody.get(field, "") if antibody else ""
            if not candidate:
                continue
            identity, coverage = alignment_metrics(amino_acids, candidate)
            if identity >= 45 and coverage >= 70:
                results.append(
                    {
                        "id": antibody_id,
                        "field": field,
                        "identity": identity,
                        "coverage": coverage,
                    }
                )
        return sorted(results, key=lambda item: (-item["identity"], -item["coverage"], item["id"]))


def alignment_metrics(reference: str, query: str) -> tuple[float, float]:
    rows, columns = len(reference) + 1, len(query) + 1
    scores = [[0] * columns for _ in range(rows)]
    trace = [[0] * columns for _ in range(rows)]
    for row in range(1, rows):
        scores[row][0], trace[row][0] = row * -2, 2
    for column in range(1, columns):
        scores[0][column], trace[0][column] = column * -2, 3
    groups = ("STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW")
    for row in range(1, rows):
        for column in range(1, columns):
            left, right = reference[row - 1], query[column - 1]
            substitution = (
                2 if left == right else 1 if any(left in g and right in g for g in groups) else -1
            )
            options = (
                scores[row - 1][column - 1] + substitution,
                scores[row - 1][column] - 2,
                scores[row][column - 1] - 2,
            )
            best = max(range(3), key=options.__getitem__)
            scores[row][column], trace[row][column] = options[best], best + 1
    row, column, matches, aligned, covered = len(reference), len(query), 0, 0, 0
    while row or column:
        direction = trace[row][column]
        if direction == 1:
            aligned += 1
            covered += 1
            matches += reference[row - 1] == query[column - 1]
            row -= 1
            column -= 1
        elif direction == 2:
            aligned += 1
            row -= 1
        else:
            aligned += 1
            column -= 1
    identity = 100 * matches / aligned if aligned else 0
    coverage = 100 * covered / len(reference) if reference else 0
    return round(identity, 3), round(coverage, 3)


def merge_rows(rows: list[dict]) -> list[dict]:
    merged: dict[str, dict] = {}
    for row in rows:
        antibody_id = row.get("antibody", {}).get("id", "")
        if antibody_id not in merged:
            merged[antibody_id] = row
            continue
        current = merged[antibody_id]
        for field in ("relationships", "evidence", "sources"):
            current[field] = sorted(set(current.get(field, [])) | set(row.get(field, [])))
        interactions = {item["id"]: item for item in current.get("interactions", [])}
        interactions.update({item["id"]: item for item in row.get("interactions", [])})
        current["interactions"] = list(interactions.values())
    return list(merged.values())


def result_payload(snapshot: Snapshot, query: dict, rows: list[dict]) -> dict:
    return {
        "schema_version": SUPPORTED_SCHEMA,
        "snapshot": snapshot.manifest.get("snapshot"),
        "query": query,
        "count": len(rows),
        "results": rows,
    }


def render(snapshot: Snapshot, payload: dict, output_format: str) -> str:
    if output_format == "json":
        return json.dumps(payload, indent=2, ensure_ascii=False)
    full = [snapshot.antibody(row["antibody"]["id"]) for row in payload["results"]]
    full = [record for record in full if record]
    if output_format == "fasta":
        entries = []
        for antibody in full:
            for chain, value in (("VH", antibody.get("heavy")), ("VL", antibody.get("light"))):
                if value:
                    entries.append(
                        f">pairs_id={antibody['id']}|chain={chain}|sequence=amino_acid\n{value}"
                    )
        return "\n".join(entries)
    handle = io.StringIO(newline="")
    fields = ["snapshot", "query", "antibody_id", "name", "sources", "heavy", "light"]
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    for antibody in full:
        writer.writerow(
            {
                "snapshot": payload.get("snapshot", ""),
                "query": payload.get("query", {}).get("value", ""),
                "antibody_id": antibody["id"],
                "name": antibody.get("name", ""),
                "sources": ";".join(antibody.get("sources", [])),
                "heavy": antibody.get("heavy", ""),
                "light": antibody.get("light", ""),
            }
        )
    return handle.getvalue()


def target_command(snapshot: Snapshot, args) -> dict:
    target = snapshot.resolve_target(args.query)
    if not target:
        raise ValueError(f"no exact target identity in this snapshot: {args.query}")
    targets = [target, *(snapshot.descendants(target) if args.include_descendants else [])]
    families = ["positive"]
    if args.include_functional:
        families.append("functional")
    if args.include_negative:
        families.append("negative")
    rows = merge_rows(
        [
            row
            for item in targets
            for family in families
            for row in snapshot.target_rows(item, family)
        ]
    )
    if args.paired:
        rows = [
            row
            for row in rows
            if row.get("antibody", {}).get("has_heavy") and row.get("antibody", {}).get("has_light")
        ]
    return result_payload(
        snapshot,
        {
            "type": "target",
            "value": args.query,
            "target_id": target["id"],
            "scope": "include_descendants" if args.include_descendants else "exact",
            "families": families,
        },
        rows,
    )


def sequence_command(snapshot: Snapshot, args) -> dict:
    first = sequence(args.sequence)
    second = sequence(args.light_sequence) if args.light_sequence else ""
    first_field = "light" if args.chain == "light" else "heavy"
    first_hits = (
        snapshot.similarity(first, first_field)
        if args.similarity
        else snapshot.exact_sequence(first, first_field)
    )
    if second:
        second_hits = (
            snapshot.similarity(second, "light")
            if args.similarity
            else snapshot.exact_sequence(second, "light")
        )
        second_ids = {hit["id"] for hit in second_hits}
        hits = [hit for hit in first_hits if hit["id"] in second_ids]
    else:
        hits = first_hits
    rows = []
    for hit in hits:
        antibody = snapshot.antibody(hit["id"])
        rows.append(
            {
                "antibody": {
                    "id": hit["id"],
                    "name": antibody.get("name", hit["id"]) if antibody else hit["id"],
                },
                "sources": antibody.get("sources", []) if antibody else [],
                "relationships": [],
                "evidence": [],
                "interactions": [],
                "similarity": hit if args.similarity else None,
            }
        )
    return result_payload(
        snapshot,
        {
            "type": "sequence",
            "value": hashlib.sha256(first.encode()).hexdigest(),
            "chain": first_field,
            "paired": bool(second),
            "similarity": args.similarity,
        },
        rows,
    )


def antibody_command(snapshot: Snapshot, args) -> dict:
    antibody = snapshot.antibody(args.id)
    if not antibody:
        raise ValueError(f"antibody ID is not present in this snapshot: {args.id}")
    return result_payload(
        snapshot,
        {"type": "antibody", "value": args.id},
        [{"antibody": {"id": args.id, "name": antibody.get("name", args.id)}}],
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(prog="pairs", description="Query a local PAIRS snapshot")
    result.add_argument("--data", type=Path, default=Path("data/v4"))
    result.add_argument("--format", choices=("json", "csv", "fasta"), default="json")
    commands = result.add_subparsers(dest="command", required=True)
    target = commands.add_parser("target")
    target.add_argument("query")
    target.add_argument("--paired", action="store_true")
    target.add_argument("--include-descendants", action="store_true")
    target.add_argument("--include-functional", action="store_true")
    target.add_argument("--include-negative", action="store_true")
    antibody = commands.add_parser("antibody")
    antibody.add_argument("id")
    sequence_parser = commands.add_parser("sequence")
    sequence_parser.add_argument("sequence")
    sequence_parser.add_argument("--chain", choices=("heavy", "light"), default="heavy")
    sequence_parser.add_argument("--light-sequence", default="")
    sequence_parser.add_argument("--similarity", action="store_true")
    return result


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    try:
        snapshot = Snapshot(args.data)
        payload = {
            "target": target_command,
            "sequence": sequence_command,
            "antibody": antibody_command,
        }[args.command](snapshot, args)
        print(render(snapshot, payload, args.format))
        return 0
    except (FileNotFoundError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
