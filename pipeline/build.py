from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
import tempfile
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Iterable

from .model import AntibodyObservation, InteractionObservation, digest, split_values, text
from .sources import ADAPTERS
from .targets import TargetResolver

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config"
SCHEMA_VERSION = 2
APP_VERSION = "2.0.0"
DATA_SUBDIR = f"v{SCHEMA_VERSION}"
USER_AGENT = "PAIRS/2.0 (Pan-Antibody Integrated Retrieval System; static scientific index)"
TARGET_PAGE_SIZE = 250


def load_sources() -> dict:
    return json.loads((CONFIG / "sources.json").read_text(encoding="utf-8"))


def _request(url: str, timeout: int = 120):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    return urllib.request.urlopen(request, timeout=timeout)


def download(url: str, dest: Path, timeout: int = 120, attempts: int = 4) -> dict:
    """Download atomically with bounded exponential retry/backoff."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    last_error: Exception | None = None

    for attempt in range(1, attempts + 1):
        tmp_path: Path | None = None
        try:
            with tempfile.NamedTemporaryFile(
                prefix=dest.name + ".", suffix=".part", dir=dest.parent, delete=False
            ) as tmp:
                tmp_path = Path(tmp.name)
                with _request(url, timeout=timeout) as response:
                    shutil.copyfileobj(response, tmp)
                    content_type = response.headers.get("Content-Type", "")
                    modified = response.headers.get("Last-Modified", "")
            tmp_path.replace(dest)
            sha = hashlib.sha256(dest.read_bytes()).hexdigest()
            return {
                "bytes": dest.stat().st_size,
                "sha256": sha,
                "content_type": content_type,
                "last_modified": modified,
                "seconds": round(time.time() - started, 2),
                "attempts": attempt,
            }
        except Exception as exc:  # pragma: no cover - network behavior varies
            last_error = exc
            if tmp_path and tmp_path.exists():
                tmp_path.unlink(missing_ok=True)
            if attempt == attempts:
                break
            time.sleep(min(2 ** (attempt - 1), 8))

    assert last_error is not None
    raise last_error


def _download_text(url: str, timeout: int = 30, attempts: int = 3) -> str:
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            with _request(url, timeout=timeout) as response:
                return response.read().decode("utf-8", errors="replace")
        except Exception as exc:  # pragma: no cover - network behavior varies
            last_error = exc
            if attempt < attempts:
                time.sleep(min(2 ** (attempt - 1), 4))
    assert last_error is not None
    raise last_error


def resolve_download_url(source_config: dict) -> tuple[str, dict]:
    """Resolve the current upstream download link from its homepage.

    Public V2 sources use discovery deliberately: a rotated upstream filename should either resolve
    to the new file or fail loudly instead of silently pinning an old snapshot forever.
    """
    fallback = source_config.get("download_url", "")
    if not source_config.get("discover_on_homepage"):
        if not fallback:
            raise ValueError(f"{source_config.get('name', 'source')} has no download URL")
        return fallback, {"discovered": False}

    html = _download_text(source_config["homepage"])
    pattern = source_config.get("download_pattern")
    if not pattern:
        raise ValueError(f"{source_config.get('name', 'source')} discovery pattern is missing")
    match = re.search(pattern, html, re.IGNORECASE)
    if match:
        return urllib.parse.urljoin(source_config["homepage"], match.group(1)), {
            "discovered": True,
            "discovered_from": source_config["homepage"],
        }

    if source_config.get("discovery_required", False):
        raise RuntimeError(
            f"{source_config.get('name', 'source')} download discovery failed; upstream markup may have changed"
        )
    if fallback:
        return fallback, {"discovered": False, "discovery_fallback": True}
    raise RuntimeError(f"{source_config.get('name', 'source')} download discovery failed")


def merge_antibody(destination: dict, observation: AntibodyObservation) -> None:
    def set_if(field: str, value) -> None:
        if value and not destination.get(field):
            destination[field] = value

    set_if("name", observation.name)
    for field in [
        "heavy",
        "light",
        "cdrh3",
        "cdrl3",
        "organism",
        "format",
        "heavy_v",
        "heavy_j",
        "light_v",
        "light_j",
        "therapeutic_status",
    ]:
        set_if(field, getattr(observation, field))

    destination.setdefault("aliases", [])
    for alias in [observation.name, *observation.aliases]:
        if alias and alias not in destination["aliases"] and alias != destination.get("name"):
            destination["aliases"].append(alias)

    destination.setdefault("structures", [])
    destination["structures"] = sorted(set(destination["structures"]) | set(observation.structures))
    destination.setdefault("sources", [])
    if observation.source not in destination["sources"]:
        destination["sources"].append(observation.source)

    destination.setdefault("source_records", [])
    destination.setdefault("source_record_count", 0)
    destination.setdefault("_source_record_keys", set())
    token = (observation.source, observation.record_id)
    if token not in destination["_source_record_keys"]:
        destination["_source_record_keys"].add(token)
        destination["source_record_count"] += 1
        if len(destination["source_records"]) < 12:
            destination["source_records"].append(
                {
                    "source": observation.source,
                    "record_id": observation.record_id,
                    "reference": observation.reference,
                    "source_url": observation.source_url,
                    "metadata": observation.metadata,
                }
            )


def normalize_target(resolver: TargetResolver, raw: str, source: str):
    raw = text(raw)
    if source == "therasabdab" and "/" in raw:
        return resolver.synonym_group(raw)
    return resolver.resolve(raw)


def expand_target_terms(raw: str, source: str) -> list[str]:
    """Fan literature co-mentions out; never convert co-occurrence into aliases."""
    raw = text(raw)
    if not raw:
        return []
    if source == "plabdab" and ";" in raw:
        return split_values(raw, separators=r";")
    return [raw]


def _antibody_shard(antibody_id: str) -> str:
    if antibody_id.startswith("ab_") and len(antibody_id) >= 5:
        return antibody_id[3:5]
    return digest(antibody_id, length=2)


def _sequence_sha(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def _review_candidates(
    targets: dict[str, dict], antibody_targets: dict[str, set[str]], limit: int = 500
) -> list[dict]:
    """Suggest possible aliases for human review without ever auto-merging them."""
    pair_overlap: dict[tuple[str, str], int] = defaultdict(int)
    for target_ids in antibody_targets.values():
        ids = sorted(target_ids)
        # Very broad records are poor synonym evidence and can explode pair counts.
        if 1 < len(ids) <= 20:
            for left, right in combinations(ids, 2):
                pair_overlap[(left, right)] += 1

    candidates = []
    for (left, right), overlap in pair_overlap.items():
        if overlap < 2:
            continue
        left_count = len(targets[left]["antibodies"])
        right_count = len(targets[right]["antibodies"])
        union = left_count + right_count - overlap
        jaccard = overlap / union if union else 0.0
        containment = overlap / min(left_count, right_count) if min(left_count, right_count) else 0.0
        if jaccard < 0.45 and containment < 0.8:
            continue
        candidates.append(
            {
                "target_a": left,
                "name_a": targets[left]["name"],
                "count_a": left_count,
                "target_b": right,
                "name_b": targets[right]["name"],
                "count_b": right_count,
                "overlap": overlap,
                "jaccard": round(jaccard, 4),
                "containment": round(containment, 4),
                "action": "REVIEW_ONLY",
            }
        )
    candidates.sort(
        key=lambda item: (-item["containment"], -item["jaccard"], -item["overlap"], item["name_a"])
    )
    return candidates[:limit]


def write_indexes(
    antibodies: dict[str, dict],
    raw_interactions: Iterable[dict],
    out: Path,
    source_counts: dict[str, dict],
) -> dict:
    resolver = TargetResolver(CONFIG / "target_aliases.json")
    targets: dict[str, dict] = {}
    by_target: dict[str, dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))
    antibody_targets: dict[str, set[str]] = defaultdict(set)
    interaction_count = 0

    # Recalculate normalized interaction counts because PLAbDab literature mention groups fan out.
    for counts in source_counts.values():
        counts["interactions"] = 0

    for interaction in raw_interactions:
        source = interaction["source"]
        for raw_target in expand_target_terms(interaction.get("target_raw", ""), source):
            target_id, target_name, aliases = normalize_target(resolver, raw_target, source)
            if not target_id:
                continue

            target = targets.setdefault(
                target_id,
                {
                    "id": target_id,
                    "name": target_name,
                    "aliases": set(),
                    "sources": set(),
                    "relationships": defaultdict(int),
                    "count": 0,
                    "antibodies": set(),
                },
            )
            if target_name and len(target_name) < len(target["name"]):
                target["name"] = target_name
            target["aliases"].update(aliases)
            target["aliases"].add(raw_target)
            target["sources"].add(source)
            target["relationships"][interaction["relationship"]] += 1
            target["count"] += 1
            target["antibodies"].add(interaction["antibody_id"])

            normalized = {
                **interaction,
                "id": "ix_"
                + digest(
                    interaction["antibody_id"],
                    target_id,
                    source,
                    interaction.get("source_record_id", ""),
                    interaction["relationship"],
                    interaction.get("epitope", ""),
                    interaction.get("note", ""),
                    raw_target,
                ),
                "target_id": target_id,
                "target_name": target_name,
                "target_raw": raw_target,
            }
            by_target[target_id][interaction["antibody_id"]].append(normalized)
            antibody_targets[interaction["antibody_id"]].add(target_id)
            interaction_count += 1
            source_counts.setdefault(source, {"records": 0, "interactions": 0})["interactions"] += 1

    for antibody_id, antibody in antibodies.items():
        all_target_ids = sorted(antibody_targets.get(antibody_id, ()))
        antibody["target_count"] = len(all_target_ids)
        antibody["targets"] = [
            {"id": target_id, "name": targets[target_id]["name"]}
            for target_id in all_target_ids[:80]
        ]
        antibody.setdefault("sources", [])
        antibody.setdefault("aliases", [])
        antibody.setdefault("structures", [])
        antibody["sources"].sort()
        antibody["aliases"].sort(key=str.casefold)
        antibody.pop("_source_record_keys", None)

    out.mkdir(parents=True, exist_ok=True)
    for subdirectory in ["targets", "antibodies", "antibody-search", "sequence-search"]:
        directory = out / subdirectory
        if directory.exists():
            shutil.rmtree(directory)
        directory.mkdir(parents=True)

    target_index = []
    target_file_bytes = 0
    largest_target_page_bytes = 0
    largest_target_pages = 0

    for target_id, target in targets.items():
        results = []
        for antibody_id, interactions in by_target[target_id].items():
            antibody = antibodies[antibody_id]
            relationships = sorted({item["relationship"] for item in interactions})
            evidence = sorted({item["evidence"] for item in interactions})
            sources = sorted({item["source"] for item in interactions})
            antibody_summary = {
                "id": antibody_id,
                "name": antibody.get("name", ""),
                "aliases": antibody.get("aliases", [])[:6],
                "organism": antibody.get("organism", ""),
                "format": antibody.get("format", ""),
                "has_heavy": bool(antibody.get("heavy")),
                "has_light": bool(antibody.get("light")),
                "heavy_length": len(antibody.get("heavy", "")),
                "light_length": len(antibody.get("light", "")),
                "structures": antibody.get("structures", [])[:12],
                "therapeutic_status": antibody.get("therapeutic_status", ""),
                "shard": _antibody_shard(antibody_id),
            }
            results.append(
                {
                    "antibody": antibody_summary,
                    "relationships": relationships,
                    "evidence": evidence,
                    "sources": sources,
                    "interactions": interactions,
                }
            )

        evidence_rank = {
            "STRUCTURE": 6,
            "MEASURED": 5,
            "CURATED": 4,
            "PUBLICATION": 3,
            "LITERATURE_METADATA": 2,
            "PATENT": 1,
            "METADATA": 0,
        }

        def result_rank(result: dict) -> int:
            evidence_score = max(
                (evidence_rank.get(item, 0) for item in result["evidence"]), default=0
            )
            return (
                evidence_score * 100
                + len(result["sources"]) * 10
                + (5 if result["antibody"].get("structures") else 0)
                + (3 if result["antibody"].get("therapeutic_status") else 0)
            )

        results.sort(
            key=lambda result: (-result_rank(result), result["antibody"].get("name", "").casefold())
        )
        target_dir_name = digest(target_id, length=20)
        target_dir = out / "targets" / target_dir_name
        target_dir.mkdir(parents=True)
        pages = []
        for page_number, start in enumerate(range(0, len(results), TARGET_PAGE_SIZE), 1):
            page_name = f"page-{page_number:03d}.json"
            page_path = target_dir / page_name
            payload = results[start : start + TARGET_PAGE_SIZE]
            page_path.write_text(
                json.dumps(payload, separators=(",", ":"), ensure_ascii=False), encoding="utf-8"
            )
            page_bytes = page_path.stat().st_size
            target_file_bytes += page_bytes
            largest_target_page_bytes = max(largest_target_page_bytes, page_bytes)
            pages.append({"file": page_name, "count": len(payload), "bytes": page_bytes})

        largest_target_pages = max(largest_target_pages, len(pages))
        target_stats = {
            "unique_results": len(results),
            "paired": sum(
                bool(result["antibody"]["has_heavy"] and result["antibody"]["has_light"])
                for result in results
            ),
            "therapeutic": sum(bool(result["antibody"]["therapeutic_status"]) for result in results),
            "structure": sum(bool(result["antibody"]["structures"]) for result in results),
            "negative": sum(
                any("does_not" in relationship or "not_" in relationship for relationship in result["relationships"])
                for result in results
            ),
        }
        target_public = {
            "id": target_id,
            "name": target["name"],
            "aliases": sorted(set(target["aliases"]), key=str.casefold),
            "count": target["count"],
            "result_count": len(results),
            "sources": sorted(target["sources"]),
            "relationships": dict(sorted(target["relationships"].items())),
            "dir": target_dir_name,
            "page_count": len(pages),
            "page_size": TARGET_PAGE_SIZE,
            "stats": target_stats,
        }
        target_index.append(target_public)
        (target_dir / "index.json").write_text(
            json.dumps({"target": target_public, "pages": pages}, separators=(",", ":"), ensure_ascii=False),
            encoding="utf-8",
        )

    target_index.sort(key=lambda item: (item["name"].casefold(), item["id"]))
    (out / "targets.json").write_text(
        json.dumps(target_index, separators=(",", ":"), ensure_ascii=False), encoding="utf-8"
    )

    antibody_search: dict[str, dict[str, dict]] = defaultdict(dict)
    antibody_shards: dict[str, dict[str, dict]] = defaultdict(dict)
    sequence_search: dict[str, dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))

    for antibody_id, antibody in antibodies.items():
        shard = _antibody_shard(antibody_id)
        search_item = {
            "id": antibody_id,
            "name": antibody.get("name", ""),
            "aliases": antibody.get("aliases", [])[:12],
            "sources": antibody.get("sources", []),
            "targets": antibody.get("targets", [])[:12],
            "paired": bool(antibody.get("heavy") and antibody.get("light")),
            "shard": shard,
        }
        terms = [antibody.get("name", ""), *antibody.get("aliases", [])]
        buckets = set()
        for term in terms:
            compact = "".join(character.lower() for character in term if character.isalnum())
            if compact:
                buckets.add((compact[:2] if len(compact) >= 2 else compact + "_")[:2])
        for bucket in buckets or {"__"}:
            antibody_search[bucket][antibody_id] = search_item
        antibody_shards[shard][antibody_id] = antibody

        for field in ["heavy", "light", "cdrh3", "cdrl3"]:
            sequence_value = antibody.get(field, "")
            if not sequence_value:
                continue
            sequence_hash = _sequence_sha(sequence_value)
            sequence_search[sequence_hash[:2]][sequence_hash].append(
                {**search_item, "field": field, "length": len(sequence_value)}
            )

    antibody_search_bytes = 0
    for bucket, payload in antibody_search.items():
        records = sorted(payload.values(), key=lambda item: (item["name"].casefold(), item["id"]))
        path = out / "antibody-search" / f"{bucket}.json"
        path.write_text(
            json.dumps(records, separators=(",", ":"), ensure_ascii=False), encoding="utf-8"
        )
        antibody_search_bytes += path.stat().st_size

    for shard, payload in antibody_shards.items():
        (out / "antibodies" / f"{shard}.json").write_text(
            json.dumps(payload, separators=(",", ":"), ensure_ascii=False), encoding="utf-8"
        )

    sequence_search_bytes = 0
    for bucket, payload in sequence_search.items():
        path = out / "sequence-search" / f"{bucket}.json"
        path.write_text(
            json.dumps(payload, separators=(",", ":"), ensure_ascii=False), encoding="utf-8"
        )
        sequence_search_bytes += path.stat().st_size

    review = _review_candidates(targets, antibody_targets)
    (out / "target-review.json").write_text(
        json.dumps(
            {
                "generated": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
                "description": "Human-review suggestions based on overlapping antibody sets. No pair is auto-merged.",
                "candidates": review,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    return {
        "antibodies": len(antibodies),
        "targets": len(targets),
        "interactions": interaction_count,
        "source_counts": dict(source_counts),
        "files": {
            "targets_index_bytes": (out / "targets.json").stat().st_size,
            "target_page_bytes": target_file_bytes,
            "largest_target_page_bytes": largest_target_page_bytes,
            "largest_target_page_count": largest_target_pages,
            "antibody_search_index_bytes": antibody_search_bytes,
            "sequence_search_index_bytes": sequence_search_bytes,
        },
        "target_review_candidates": len(review),
    }


def compile_data(
    source_paths: dict[str, Path],
    source_config: dict,
    out: Path,
    max_records: int | None = None,
) -> dict:
    antibodies: dict[str, dict] = {}
    raw_interactions: list[dict] = []
    source_counts: dict[str, dict] = defaultdict(lambda: {"records": 0, "interactions": 0})

    for source, path in source_paths.items():
        config = source_config[source]
        adapter = ADAPTERS[config["adapter"]]
        for index, (observation, interactions) in enumerate(adapter(path, config.get("homepage", ""))):
            if max_records is not None and index >= max_records:
                break
            antibody_id = observation.identity()
            antibody = antibodies.setdefault(antibody_id, {"id": antibody_id, "name": observation.name})
            merge_antibody(antibody, observation)
            source_counts[source]["records"] += 1
            for interaction in interactions:
                raw_interactions.append(
                    {
                        "antibody_id": antibody_id,
                        "target_raw": interaction.target_raw,
                        "relationship": interaction.relationship,
                        "evidence": interaction.evidence,
                        "source": interaction.source,
                        "source_record_id": interaction.source_record_id,
                        "reference": interaction.reference,
                        "epitope": interaction.epitope,
                        "assay": interaction.assay,
                        "note": interaction.note,
                    }
                )

    return write_indexes(antibodies, raw_interactions, out, source_counts)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Build the static PAIRS dataset")
    parser.add_argument("--output", default=str(ROOT / "data" / DATA_SUBDIR))
    parser.add_argument("--cache", default=str(ROOT / ".cache" / "sources"))
    parser.add_argument("--offline", action="store_true", help="Use source files already present in --cache")
    parser.add_argument(
        "--max-records", type=int, default=None, help="Limit records per source (for local/demo builds)"
    )
    parser.add_argument(
        "--source", action="append", dest="sources", help="Build only selected source key; repeatable"
    )
    parser.add_argument(
        "--allow-partial", action="store_true", help="Continue when a source download fails"
    )
    args = parser.parse_args(argv)

    config = load_sources()
    selected = args.sources or [key for key, value in config.items() if value.get("enabled_public")]
    cache = Path(args.cache)
    cache.mkdir(parents=True, exist_ok=True)
    source_paths: dict[str, Path] = {}
    statuses: dict[str, dict] = {}

    for source in selected:
        if source not in config:
            raise SystemExit(f"Unknown source: {source}")
        source_config = config[source]
        if source_config.get("adapter") not in ADAPTERS:
            print(f"skip {source}: no enabled adapter", file=sys.stderr)
            continue
        path = cache / f"{source}.csv"
        try:
            info: dict = {}
            if not args.offline:
                print(f"Downloading {source} ...", flush=True)
                resolved_url, discovery = resolve_download_url(source_config)
                info = download(resolved_url, path)
                info.update(discovery)
                info["download_url"] = resolved_url
            elif not path.exists():
                raise FileNotFoundError(path)
            source_paths[source] = path
            statuses[source] = {
                "ok": True,
                "name": source_config["name"],
                "license": source_config.get("license", ""),
                "homepage": source_config.get("homepage", ""),
                **info,
            }
        except Exception as exc:
            statuses[source] = {
                "ok": False,
                "name": source_config["name"],
                "error": str(exc),
                "homepage": source_config.get("homepage", ""),
            }
            print(f"WARNING {source}: {exc}", file=sys.stderr)
            if not args.allow_partial:
                raise

    if not source_paths:
        raise SystemExit("No sources available to build")

    out = Path(args.output)
    stats = compile_data(source_paths, config, out, args.max_records)
    for source, counts in stats.get("source_counts", {}).items():
        statuses.setdefault(source, {}).update(counts)

    snapshot = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "app_version": APP_VERSION,
        "data_path": f"data/{DATA_SUBDIR}",
        "snapshot": snapshot,
        "stats": stats,
        "sources": statuses,
        "sources_expected": len(selected),
        "sources_ok": sum(1 for status in statuses.values() if status.get("ok")),
        "disclaimer": (
            "Public scientific metadata aggregator. Verify sequence identity, target assignment and "
            "experimental evidence against the cited primary source before experimental use."
        ),
    }
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(manifest["stats"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
