from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import shutil
import sys
import tempfile
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Iterable

from .model import AntibodyObservation, digest, split_values, text
from .search import (
    build_exact_antibody_index,
    normalize_search_term,
    search_bucket,
)
from .sources import ADAPTERS
from .targets import TargetResolver

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config"
SCHEMA_VERSION = 4
DATA_CONTRACT_REVISION = 8
APP_VERSION = "4.0.0"
DATA_SUBDIR = f"v{SCHEMA_VERSION}"
USER_AGENT = "PAIRS/4.0 (Pan-Antibody Integrated Retrieval System; static scientific index)"
TARGET_PAGE_SIZE = 100
SIMILARITY_K = 5
SIMILARITY_SIGNATURE_SIZE = 32
SIMILARITY_MAX_POSTING = 250
CLUSTER_THRESHOLDS = (99, 95, 90)
CLUSTER_MIN_COVERAGE = 90
CLUSTER_MIN_SHARED_SIGNATURES = 5

# Keep the retrieval contract explicit.  A target page is a claim about a
# sequence, so literature co-occurrence and negative observations must never
# become default positive hits. Functional observations remain searchable as
# positive evidence, but are retained in their own page family/collection.
DIRECT_POSITIVE_RELATIONSHIPS = frozenset({"binds", "targets"})
FUNCTIONAL_POSITIVE_RELATIONSHIPS = frozenset({"neutralizes", "protects"})
NEGATIVE_RELATIONSHIP_PREFIXES = ("does_not_", "not_")
LITERATURE_RELATIONSHIPS = frozenset({"mentioned_with"})


def _is_negative_relationship(relationship: str) -> bool:
    return relationship.startswith(NEGATIVE_RELATIONSHIP_PREFIXES)


def load_sources() -> dict:
    return json.loads((CONFIG / "sources.json").read_text(encoding="utf-8"))


def _request(url: str, timeout: int = 120, headers: dict | None = None):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, **(headers or {})})
    return urllib.request.urlopen(request, timeout=timeout)


def download(
    url: str,
    dest: Path,
    timeout: int = 120,
    attempts: int = 4,
    accept: str = "",
) -> dict:
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
                with _request(
                    url,
                    timeout=timeout,
                    headers={"Accept": accept} if accept else None,
                ) as response:
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


def download_postgrest_csv(
    url: str,
    dest: Path,
    timeout: int = 120,
    attempts: int = 4,
    page_size: int = 10_000,
    order: str = "receptor_group_id.asc,epitope__iedb_iri.asc,assay__iedb_ids.asc",
) -> dict:
    """Download a complete PostgREST CSV endpoint using deterministic pagination."""
    started = time.time()
    dest.parent.mkdir(parents=True, exist_ok=True)
    temporary_handle = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        dir=dest.parent,
        delete=False,
    )
    temporary = Path(temporary_handle.name)
    temporary_handle.close()
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            expected_header = None
            offset = 0
            page_count = 0
            while True:
                separator = "&" if "?" in url else "?"
                page_url = f"{url}{separator}limit={page_size}&offset={offset}&order={order}"
                last_error = None
                for attempt in range(1, attempts + 1):
                    try:
                        with _request(
                            page_url,
                            timeout=timeout,
                            headers={"Accept": "text/csv"},
                        ) as response:
                            payload = response.read().decode("utf-8-sig")
                        break
                    except Exception as exc:  # pragma: no cover - network behavior varies
                        last_error = exc
                        if attempt == attempts:
                            raise
                        time.sleep(min(2 ** (attempt - 1), 8))
                else:  # pragma: no cover - defensive; loop either breaks or raises
                    raise RuntimeError(f"download failed for {page_url}: {last_error}")
                rows = list(csv.reader(io.StringIO(payload, newline="")))
                if not rows:
                    break
                header, data_rows = rows[0], rows[1:]
                if expected_header is None:
                    expected_header = header
                    writer.writerow(header)
                elif header != expected_header:
                    raise ValueError("PostgREST CSV header changed between pages")
                writer.writerows(data_rows)
                page_count += 1
                offset += len(data_rows)
                if len(data_rows) < page_size:
                    break
                if page_count >= 100:
                    raise RuntimeError("PostgREST pagination exceeded 100 pages")
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    temporary.replace(dest)
    payload = dest.read_bytes()
    return {
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "content_type": "text/csv; charset=utf-8",
        "seconds": round(time.time() - started, 2),
        "attempts": attempts,
        "pages": page_count,
        "rows": offset,
    }


def _postgres_array(value: str) -> list[str]:
    value = (value or "").strip()
    if not value or value == "{}":
        return []
    return [item.strip().strip('"') for item in value.strip("{}").split(",") if item.strip()]


def download_iedb_support(cache: Path) -> dict:
    search_path = cache / "iedb-bcr-search.csv"
    search_info = download_postgrest_csv(
        "https://query-api.iedb.org/bcr_search",
        search_path,
        order="receptor_group_id.asc",
    )
    search_rows = list(_open_csv_for_build(search_path))
    if search_rows and not {"receptor_group_id", "iedb_assay_ids"}.issubset(search_rows[0]):
        raise ValueError("IEDB bcr_search support schema is missing linkage columns")
    assay_ids = {
        assay_id
        for row in search_rows
        for assay_id in _postgres_array(row.get("iedb_assay_ids", ""))
    }
    invalid_ids = sorted(assay_id for assay_id in assay_ids if not assay_id.isdigit())
    if invalid_ids:
        raise ValueError(f"IEDB returned nonnumeric assay IDs: {invalid_ids[:3]}")
    assay_ids = sorted(assay_ids, key=int)
    measurement_path = cache / "iedb-bcell-export.csv"
    temporary = measurement_path.with_suffix(".csv.part")
    header = None
    row_count = 0

    def fetch_assay_batch(batch: list[str]) -> list[list[str]]:
        url = (
            "https://query-api.iedb.org/bcell_export?assay_id=in.("
            + ",".join(batch)
            + ")&order=assay_id.asc"
        )
        last_error = None
        for attempt in range(1, 5):
            try:
                with _request(url, timeout=60, headers={"Accept": "text/csv"}) as response:
                    return list(csv.reader(io.StringIO(response.read().decode("utf-8-sig"))))
            except Exception as exc:  # pragma: no cover - upstream network behavior
                last_error = exc
                if attempt < 4:
                    time.sleep(min(2 ** (attempt - 1), 8))
        assert last_error is not None
        raise last_error

    batches = [assay_ids[start : start + 100] for start in range(0, len(assay_ids), 100)]
    try:
        with temporary.open("w", encoding="utf-8", newline="") as output:
            writer = csv.writer(output)
            with ThreadPoolExecutor(max_workers=4) as executor:
                batch_rows = executor.map(fetch_assay_batch, batches)
                for rows in batch_rows:
                    if not rows:
                        continue
                    if header is None:
                        header = rows[0]
                        required = {"assay_id", "assay__method"}
                        if not required.issubset(header):
                            raise ValueError(
                                "IEDB bcell_export support schema is missing required columns"
                            )
                        writer.writerow(header)
                    elif rows[0] != header:
                        raise ValueError("IEDB bcell_export schema changed between assay batches")
                    writer.writerows(rows[1:])
                    row_count += len(rows) - 1
        if assay_ids and not row_count:
            raise ValueError("IEDB returned no assay rows for linked receptor assays")
        temporary.replace(measurement_path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return {
        "bcr_search": search_info,
        "linked_assay_ids": len(assay_ids),
        "bcell_export_rows": row_count,
    }


def inspect_cached_iedb_support(cache: Path) -> dict:
    """Describe cached IEDB linkage files without contacting the upstream API."""

    search_path = cache / "iedb-bcr-search.csv"
    search_rows = list(_open_csv_for_build(search_path))
    assay_ids = {
        assay_id
        for row in search_rows
        for assay_id in _postgres_array(row.get("iedb_assay_ids", ""))
    }
    measurement_path = cache / "iedb-bcell-export.csv"
    with measurement_path.open(encoding="utf-8-sig", newline="") as handle:
        measurement_rows = max(0, sum(1 for _ in csv.DictReader(handle)))
    return {
        "bcr_search": {
            "bytes": search_path.stat().st_size,
            "sha256": hashlib.sha256(search_path.read_bytes()).hexdigest(),
            "content_type": "text/csv; charset=utf-8",
            "rows": len(search_rows),
        },
        "linked_assay_ids": len(assay_ids),
        "bcell_export_rows": measurement_rows,
    }


def _open_csv_for_build(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        yield from csv.DictReader(handle)


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
    for field in ["vh_nt_source", "vl_nt_source"]:
        incoming = getattr(observation, field)
        existing = destination.get(field, "")
        if incoming and existing and incoming != existing:
            conflict = {
                "field": field,
                "existing_value": existing,
                "incoming_value": incoming,
                "incoming_source": observation.source,
                "incoming_source_record_id": observation.record_id,
            }
            destination.setdefault("source_conflicts", [])
            if conflict not in destination["source_conflicts"]:
                destination["source_conflicts"].append(conflict)
    destination.setdefault("metadata", {})
    for key, value in observation.metadata.items():
        if value and key not in destination["metadata"]:
            destination["metadata"][key] = value
    for field in [
        "heavy",
        "light",
        "vh_nt_source",
        "vl_nt_source",
        "cdrh3",
        "cdrl3",
        "cdrh1",
        "cdrh2",
        "cdrl1",
        "cdrl2",
        "organism",
        "format",
        "heavy_v",
        "heavy_d",
        "heavy_j",
        "light_v",
        "light_j",
        "therapeutic_status",
    ]:
        set_if(field, getattr(observation, field))

    if observation.nucleotide_provenance:
        destination.setdefault("nucleotide_provenance", {})
        for chain, provenance in observation.nucleotide_provenance.items():
            destination["nucleotide_provenance"].setdefault(chain, provenance)
    destination.setdefault("source_nucleotide_records", [])
    nucleotide_keys = {
        (
            record.get("source", ""),
            record.get("source_record_id", ""),
            record.get("chain", ""),
            record.get("sequence", ""),
            record.get("scope", "unknown"),
        )
        for record in destination["source_nucleotide_records"]
    }
    for record in observation.source_nucleotide_records:
        normalized_record = {
            "source": record.get("source", observation.source),
            "source_record_id": record.get("source_record_id", observation.record_id),
            "chain": record.get("chain", ""),
            "sequence": record.get("sequence", ""),
            "scope": record.get("scope", "unknown"),
            "source_field": record.get("source_field", ""),
        }
        key = (
            normalized_record["source"],
            normalized_record["source_record_id"],
            normalized_record["chain"],
            normalized_record["sequence"],
            normalized_record["scope"],
        )
        if key not in nucleotide_keys:
            destination["source_nucleotide_records"].append(normalized_record)
            nucleotide_keys.add(key)
    if observation.chain_annotations:
        destination.setdefault("chain_annotations", {})
        for chain, annotation in observation.chain_annotations.items():
            if chain not in destination["chain_annotations"]:
                destination["chain_annotations"][chain] = annotation
            elif destination["chain_annotations"][chain] != annotation:
                destination.setdefault("source_conflicts", []).append(
                    {
                        "field": f"chain_annotations.{chain}",
                        "existing_value": destination["chain_annotations"][chain],
                        "incoming_value": annotation,
                        "incoming_source": observation.source,
                        "incoming_source_record_id": observation.record_id,
                    }
                )

    destination.setdefault("aliases", [])
    for alias in [observation.name, *observation.aliases]:
        if alias and alias not in destination["aliases"] and alias != destination.get("name"):
            destination["aliases"].append(alias)

    destination.setdefault("structures", [])
    destination["structures"] = sorted(set(destination["structures"]) | set(observation.structures))
    destination.setdefault("structure_tiers", {})
    for tier, structures in observation.structure_tiers.items():
        destination["structure_tiers"][tier] = sorted(
            set(destination["structure_tiers"].get(tier, [])) | set(structures)
        )
    if observation.construct:
        destination.setdefault("constructs", [])
        construct_id = observation.construct.get("id")
        if construct_id and not any(
            item.get("id") == construct_id for item in destination["constructs"]
        ):
            destination["constructs"].append(observation.construct)
    if observation.arm:
        destination.setdefault("arms", [])
        arm_id = observation.arm.get("id")
        if arm_id and not any(item.get("id") == arm_id for item in destination["arms"]):
            destination["arms"].append(observation.arm)
    destination.setdefault("sources", [])
    if observation.source not in destination["sources"]:
        destination["sources"].append(observation.source)

    destination.setdefault("source_records", [])
    destination.setdefault("source_record_count", 0)
    destination.setdefault("_source_record_keys", set())
    token = (observation.source, observation.record_id)
    existing_source_record = next(
        (
            record
            for record in destination["source_records"]
            if (record.get("source"), record.get("record_id")) == token
        ),
        None,
    )
    if existing_source_record is not None and observation.source_nucleotide_records:
        existing_source_record.setdefault("nucleotide_records", [])
        existing_keys = {
            (
                record.get("source", ""),
                record.get("source_record_id", ""),
                record.get("chain", ""),
                record.get("sequence", ""),
                record.get("scope", "unknown"),
            )
            for record in existing_source_record["nucleotide_records"]
        }
        for record in observation.source_nucleotide_records:
            key = (
                record.get("source", ""),
                record.get("source_record_id", ""),
                record.get("chain", ""),
                record.get("sequence", ""),
                record.get("scope", "unknown"),
            )
            if key not in existing_keys:
                existing_source_record["nucleotide_records"].append(dict(record))
                existing_keys.add(key)
    if token not in destination["_source_record_keys"]:
        destination["_source_record_keys"].add(token)
        destination["source_record_count"] += 1
        record_date_field = next(
            (
                field
                for field in ("date_added", "added", "update_date")
                if observation.metadata.get(field)
            ),
            "",
        )
        destination["source_records"].append(
            {
                "source": observation.source,
                "record_id": observation.record_id,
                "reference": observation.reference,
                "source_url": observation.source_url,
                "record_url": observation.record_url,
                "link_scope": observation.link_scope,
                "metadata": observation.metadata,
                "record_date": observation.metadata.get(record_date_field, ""),
                "record_date_field": record_date_field,
                "nucleotide_records": [
                    dict(record) for record in observation.source_nucleotide_records
                ],
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


def _fnv1a_32(value: str) -> int:
    result = 2166136261
    for character in value:
        result ^= ord(character)
        result = (result * 16777619) & 0xFFFFFFFF
    return result


def _sequence_signature(sequence: str, k: int = SIMILARITY_K) -> list[str]:
    """Bottom-k hashed peptide signature used only for candidate retrieval."""
    if len(sequence) < k:
        return []
    hashes = {
        _fnv1a_32(sequence[index : index + k])
        for index in range(len(sequence) - k + 1)
        if set(sequence[index : index + k]) <= set("ACDEFGHIKLMNPQRSTVWY")
    }
    return [f"{value:08x}" for value in sorted(hashes)[:SIMILARITY_SIGNATURE_SIZE]]


def sequence_quality(antibody: dict) -> dict:
    heavy = antibody.get("heavy", "")
    light = antibody.get("light", "")
    pairing = (
        "paired"
        if heavy and light
        else "heavy_only"
        if heavy
        else "light_only"
        if light
        else "no_full_chain"
    )
    ambiguous = sorted(set(heavy + light) & set("*BJOUXZ"))
    metadata = antibody.get("metadata", {})
    descriptor = f"{antibody.get('format', '')} {metadata.get('domain_type', '')}".casefold()
    explicit_vhh = bool(
        heavy
        and not light
        and re.search(
            r"(?:\bvhh\b|\bnb\b|nanob(?:ody|odies)|single[- ]domain|single[- ]chain\s+vhh)",
            descriptor,
        )
    )
    return {
        "pairing": pairing,
        "heavy_length": len(heavy),
        "light_length": len(light),
        "explicit_vhh": explicit_vhh,
        "ambiguous_residues": ambiguous,
        "source_format_quarantined": bool(metadata.get("sequence_quarantine")),
        "completeness": metadata.get("sequence_completeness") or "unknown_not_inferred",
        "source_nucleotide_available": bool(
            antibody.get("vh_nt_source") or antibody.get("vl_nt_source")
        ),
    }


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
        # Alias suggestions must be based on positive target assignments, not
        # negative observations or literature context.
        left_count = len(targets[left].get("direct_antibodies", set()))
        right_count = len(targets[right].get("direct_antibodies", set()))
        union = left_count + right_count - overlap
        jaccard = overlap / union if union else 0.0
        containment = (
            overlap / min(left_count, right_count) if min(left_count, right_count) else 0.0
        )
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


def write_similarity_indexes(antibodies: dict[str, dict], out: Path) -> dict:
    root = out / "similarity"
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    indexes = {
        "heavy": defaultdict(lambda: defaultdict(set)),
        "light": defaultdict(lambda: defaultdict(set)),
    }
    for antibody_id, antibody in antibodies.items():
        for field in indexes:
            sequence_value = antibody.get(field, "")
            for signature_hash in _sequence_signature(sequence_value):
                indexes[field][signature_hash[0]][signature_hash].add(antibody_id)

    stats = {}
    for field, buckets in indexes.items():
        directory = root / field
        directory.mkdir(parents=True)
        indexed_hashes = skipped_hashes = postings = total_bytes = 0
        for bucket in "0123456789abcdef":
            values = buckets.get(bucket, {})
            payload = {}
            for signature_hash, antibody_ids in sorted(values.items()):
                if len(antibody_ids) > SIMILARITY_MAX_POSTING:
                    skipped_hashes += 1
                    continue
                payload[signature_hash] = sorted(antibody_ids)
                indexed_hashes += 1
                postings += len(antibody_ids)
            path = directory / f"{bucket}.json"
            path.write_text(
                json.dumps(payload, separators=(",", ":"), ensure_ascii=False),
                encoding="utf-8",
            )
            total_bytes += path.stat().st_size
        metadata = {
            "version": 1,
            "field": field,
            "algorithm": "bottom-k distinct amino-acid k-mer FNV-1a signature",
            "k": SIMILARITY_K,
            "signature_size": SIMILARITY_SIGNATURE_SIZE,
            "max_posting_frequency": SIMILARITY_MAX_POSTING,
            "bucket_prefix_length": 1,
            "indexed_hashes": indexed_hashes,
            "skipped_high_frequency_hashes": skipped_hashes,
            "postings": postings,
        }
        index_path = directory / "index.json"
        index_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        total_bytes += index_path.stat().st_size
        stats[field] = {**metadata, "bytes": total_bytes}
    return stats


class _DisjointSet:
    def __init__(self, values: Iterable[str]):
        self.parent = {value: value for value in values}

    def find(self, value: str) -> str:
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left: str, right: str) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root == right_root:
            return
        first, second = sorted((left_root, right_root))
        self.parent[second] = first


def _global_edit_metrics(left: str, right: str, minimum_identity: int = 90):
    """Return deterministic global edit identity/coverage, bounded at the requested identity."""
    if not left or not right:
        return None
    longer, shorter = max(len(left), len(right)), min(len(left), len(right))
    coverage = 100 * shorter / longer
    if coverage + 1e-9 < CLUSTER_MIN_COVERAGE:
        return None
    maximum_edits = int(longer * (100 - minimum_identity) / 100 + 1e-9)
    if longer - shorter > maximum_edits:
        return None
    infinity = maximum_edits + longer + 1
    previous = list(range(len(right) + 1))
    for row, left_residue in enumerate(left, 1):
        current = [infinity] * (len(right) + 1)
        if row <= maximum_edits:
            current[0] = row
        start = max(1, row - maximum_edits)
        end = min(len(right), row + maximum_edits)
        row_minimum = infinity
        for column in range(start, end + 1):
            current[column] = min(
                previous[column] + 1,
                current[column - 1] + 1,
                previous[column - 1] + (left_residue != right[column - 1]),
            )
            row_minimum = min(row_minimum, current[column])
        if row_minimum > maximum_edits:
            return None
        previous = current
    distance = previous[len(right)]
    if distance > maximum_edits:
        return None
    return round(100 * (longer - distance) / longer, 4), round(coverage, 4)


def _cluster_scope(antibodies: dict[str, dict], scope: str) -> tuple[dict, dict]:
    fields = ("heavy", "light") if scope == "paired" else (scope,)
    records = {
        antibody_id: antibody
        for antibody_id, antibody in antibodies.items()
        if all(len(antibody.get(field, "")) >= 40 for field in fields)
    }
    signatures = {
        field: {
            antibody_id: set(_sequence_signature(antibody[field]))
            for antibody_id, antibody in records.items()
        }
        for field in fields
    }
    postings = {field: defaultdict(list) for field in fields}
    for field in fields:
        for antibody_id, values in signatures[field].items():
            for signature_hash in values:
                postings[field][signature_hash].append(antibody_id)

    sets = {threshold: _DisjointSet(records) for threshold in CLUSTER_THRESHOLDS}
    exact_groups = defaultdict(list)
    for antibody_id, antibody in records.items():
        exact_groups[tuple(antibody[field] for field in fields)].append(antibody_id)
    for members in exact_groups.values():
        for other in members[1:]:
            for disjoint in sets.values():
                disjoint.union(members[0], other)

    compared = accepted = 0
    anchor = fields[0]
    for antibody_id in sorted(records):
        counts = defaultdict(int)
        for signature_hash in signatures[anchor][antibody_id]:
            candidates = postings[anchor][signature_hash]
            if len(candidates) <= SIMILARITY_MAX_POSTING:
                for candidate_id in candidates:
                    if candidate_id > antibody_id:
                        counts[candidate_id] += 1
        for candidate_id, shared in counts.items():
            if shared < CLUSTER_MIN_SHARED_SIGNATURES:
                continue
            if (
                len(fields) == 2
                and len(signatures[fields[1]][antibody_id] & signatures[fields[1]][candidate_id])
                < CLUSTER_MIN_SHARED_SIGNATURES
            ):
                continue
            metrics = [
                _global_edit_metrics(records[antibody_id][field], records[candidate_id][field])
                for field in fields
            ]
            compared += 1
            if any(metric is None for metric in metrics):
                continue
            score = min(metric[0] for metric in metrics if metric)
            accepted += 1
            for threshold, disjoint in sets.items():
                if score + 1e-9 >= threshold:
                    disjoint.union(antibody_id, candidate_id)

    output = {}
    summary = {"records": len(records), "compared_pairs": compared, "accepted_edges": accepted}
    for threshold, disjoint in sets.items():
        groups = defaultdict(list)
        for antibody_id in sorted(records):
            groups[disjoint.find(antibody_id)].append(antibody_id)
        clusters = []
        lookup = {}
        for members in sorted((members for members in groups.values() if len(members) > 1)):
            representative = min(members)
            cluster_id = f"seq_{scope}_{threshold}_{representative[3:]}"
            clusters.append(
                {
                    "id": cluster_id,
                    "representative_id": representative,
                    "size": len(members),
                    "members": members,
                }
            )
            lookup.update({member: cluster_id for member in members})
        output[threshold] = {"clusters": clusters, "lookup": lookup}
        summary[str(threshold)] = {"clusters": len(clusters), "members": len(lookup)}
    return output, summary


def write_cluster_indexes(antibodies: dict[str, dict], out: Path) -> dict:
    root = out / "clusters"
    if root.exists():
        shutil.rmtree(root)
    metadata = {
        "version": 1,
        "label": "sequence clusters; not clonal lineages",
        "algorithm": "single-link clusters over global Levenshtein identity after bottom-k k-mer candidate retrieval",
        "identity_definition": "100 * (longer_length - global_edit_distance) / longer_length",
        "minimum_coverage": CLUSTER_MIN_COVERAGE,
        "thresholds": list(CLUSTER_THRESHOLDS),
        "candidate_min_shared_signatures": CLUSTER_MIN_SHARED_SIGNATURES,
        "scopes": {},
    }
    for scope in ("heavy", "light", "paired"):
        output, summary = _cluster_scope(antibodies, scope)
        metadata["scopes"][scope] = summary
        for threshold, payload in output.items():
            directory = root / scope / str(threshold)
            (directory / "lookup").mkdir(parents=True, exist_ok=True)
            (directory / "index.json").write_text(
                json.dumps(
                    {
                        "scope": scope,
                        "threshold": threshold,
                        "minimum_coverage": CLUSTER_MIN_COVERAGE,
                        "label": "sequence clusters; not clonal lineages",
                        "clusters": payload["clusters"],
                    },
                    separators=(",", ":"),
                ),
                encoding="utf-8",
            )
            by_shard = defaultdict(dict)
            for antibody_id, cluster_id in payload["lookup"].items():
                by_shard[antibody_id[3:5]][antibody_id] = cluster_id
            for shard in (
                left + right for left in "0123456789abcdef" for right in "0123456789abcdef"
            ):
                (directory / "lookup" / f"{shard}.json").write_text(
                    json.dumps(by_shard.get(shard, {}), separators=(",", ":")),
                    encoding="utf-8",
                )
    root.mkdir(parents=True, exist_ok=True)
    (root / "index.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata


def write_indexes(
    antibodies: dict[str, dict],
    raw_interactions: Iterable[dict],
    out: Path,
    source_counts: dict[str, dict],
    validate_examples: bool = False,
) -> dict:
    resolver = TargetResolver(CONFIG / "target_aliases.json")
    targets: dict[str, dict] = {}
    by_target: dict[str, dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))
    # These collections intentionally have different meanings.  In
    # particular, only direct positive assignments feed the alias-review set
    # and the public direct_targets list.
    antibody_direct_targets: dict[str, set[str]] = defaultdict(set)
    antibody_functional_targets: dict[str, set[str]] = defaultdict(set)
    antibody_negative_evidence: dict[str, list[dict]] = defaultdict(list)
    antibody_literature_mentions: dict[str, list[dict]] = defaultdict(list)
    interaction_count = 0
    indexed_interaction_count = 0
    literature_mention_count = 0

    # Recalculate normalized interaction counts because PLAbDab literature mention groups fan out.
    for counts in source_counts.values():
        counts["interactions"] = 0

    for interaction in raw_interactions:
        source = interaction["source"]
        for raw_target in expand_target_terms(interaction.get("target_raw", ""), source):
            target_id, target_name, aliases = normalize_target(resolver, raw_target, source)
            if not target_id:
                continue

            relationship = interaction["relationship"]
            interaction_count += 1

            # Literature context is retained on the antibody record for
            # provenance/discovery, but can never create a target, target
            # count, target page result, or direct target assignment.
            if relationship in LITERATURE_RELATIONSHIPS:
                antibody_literature_mentions[interaction["antibody_id"]].append(
                    {
                        **interaction,
                        "target_id": target_id,
                        "target_name": target_name,
                        "target_raw": raw_target,
                    }
                )
                literature_mention_count += 1
                source_counts.setdefault(source, {"records": 0, "interactions": 0})[
                    "interactions"
                ] += 1
                continue

            existing_target = targets.get(target_id)
            if (
                existing_target
                and existing_target["name"].strip().casefold() != target_name.strip().casefold()
            ):
                raise ValueError(
                    f"target ID collision: {target_id} maps to both "
                    f"{existing_target['name']!r} and {target_name!r}"
                )
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
                    "external_ids": set(),
                    "entity": resolver.entity(target_name),
                    "children": set(),
                },
            )
            if target_name and len(target_name) < len(target["name"]):
                target["name"] = target_name
            target["aliases"].update(aliases)
            target["aliases"].add(raw_target)
            target["sources"].add(source)
            if interaction.get("target_external_id"):
                target["external_ids"].add(interaction["target_external_id"])
            target["relationships"][interaction["relationship"]] += 1
            target["count"] += 1
            if relationship in DIRECT_POSITIVE_RELATIONSHIPS:
                target["direct_count"] = target.get("direct_count", 0) + 1
                target["direct_antibodies"] = target.get("direct_antibodies", set())
                target["direct_antibodies"].add(interaction["antibody_id"])
            elif relationship in FUNCTIONAL_POSITIVE_RELATIONSHIPS:
                target["functional_interaction_count"] = (
                    target.get("functional_interaction_count", 0) + 1
                )
            elif _is_negative_relationship(relationship):
                target["negative_count"] = target.get("negative_count", 0) + 1
                target["negative_antibodies"] = target.get("negative_antibodies", set())
                target["negative_antibodies"].add(interaction["antibody_id"])
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
                    interaction.get("assay", ""),
                    interaction.get("target_external_id", ""),
                    ",".join(sorted(interaction.get("assay_ids", []))),
                    interaction.get("receptor_group_id", ""),
                    json.dumps(interaction.get("measurements", []), sort_keys=True),
                    interaction.get("note", ""),
                    raw_target,
                ),
                "target_id": target_id,
                "target_name": target_name,
                "target_raw": raw_target,
            }
            by_target[target_id][interaction["antibody_id"]].append(normalized)
            indexed_interaction_count += 1
            if relationship in DIRECT_POSITIVE_RELATIONSHIPS:
                antibody_direct_targets[interaction["antibody_id"]].add(target_id)
            elif relationship in FUNCTIONAL_POSITIVE_RELATIONSHIPS:
                antibody_functional_targets[interaction["antibody_id"]].add(target_id)
            elif _is_negative_relationship(relationship):
                antibody_negative_evidence[interaction["antibody_id"]].append(normalized)
            source_counts.setdefault(source, {"records": 0, "interactions": 0})["interactions"] += 1

    for target_id, target in targets.items():
        entity = target.get("entity", {})
        parent_name = entity.get("parent")
        if not parent_name:
            continue
        parent_id, _, _ = resolver.resolve(parent_name)
        if parent_id not in targets:
            continue
        target["hierarchy"] = {
            "parent_id": parent_id,
            "relation": entity.get("relation", "related_to"),
            "scope": "curated_exact_entity",
            "provenance": entity.get("mapping_provenance", {}),
        }
        targets[parent_id]["children"].add(target_id)

    for antibody_id, antibody in antibodies.items():
        antibody["sequence_quality"] = sequence_quality(antibody)
        direct_target_ids = sorted(antibody_direct_targets.get(antibody_id, ()))
        functional_target_ids = sorted(antibody_functional_targets.get(antibody_id, ()))
        antibody["target_count"] = len(direct_target_ids)
        antibody["direct_target_count"] = len(direct_target_ids)
        antibody["functional_target_count"] = len(functional_target_ids)
        antibody["direct_targets"] = [
            {"id": target_id, "name": targets[target_id]["name"]}
            for target_id in direct_target_ids[:80]
        ]
        antibody["functional_targets"] = [
            {"id": target_id, "name": targets[target_id]["name"]}
            for target_id in functional_target_ids[:80]
        ]
        antibody["negative_evidence"] = antibody_negative_evidence.get(antibody_id, [])
        antibody["literature_mentions"] = antibody_literature_mentions.get(antibody_id, [])
        antibody.setdefault("sources", [])
        antibody.setdefault("aliases", [])
        antibody.setdefault("structures", [])
        # In the current contract ``structures`` contains only exact matches;
        # expose the meaning explicitly so the frontend never has to infer it.
        antibody["exact_structures"] = list(antibody["structures"])
        antibody["sources"].sort()
        antibody["aliases"].sort(key=str.casefold)
        antibody.pop("_source_record_keys", None)

    out.mkdir(parents=True, exist_ok=True)
    for subdirectory in [
        "targets",
        "antibodies",
        "antibody-search",
        "antibody-exact",
        "sequence-search",
        "sequence",
        "similarity",
    ]:
        directory = out / subdirectory
        if directory.exists():
            shutil.rmtree(directory)
        directory.mkdir(parents=True)

    target_index = []
    target_file_bytes = 0
    largest_target_page_bytes = 0
    largest_target_pages = 0

    for target_id, target in targets.items():

        def make_result(antibody_id: str, interactions: list[dict]) -> dict:
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
                "structure_tiers": {
                    tier: values[:12]
                    for tier, values in antibody.get("structure_tiers", {}).items()
                },
                "exact_structures": antibody.get("structures", [])[:12],
                "homologous_structures": sorted(
                    {
                        structure
                        for tier, values in antibody.get("structure_tiers", {}).items()
                        if tier != "100%" and tier not in {"exact", "exact_100"}
                        for structure in values
                    }
                )[:12],
                "therapeutic_status": antibody.get("therapeutic_status", ""),
                "sequence_quality": antibody.get("sequence_quality", {}),
                "shard": _antibody_shard(antibody_id),
            }
            return {
                "antibody": antibody_summary,
                "relationships": relationships,
                "evidence": evidence,
                "sources": sources,
                "interactions": interactions,
            }

        # The ordinary target pages are deliberately direct-positive-only.
        # Functional and negative observations are available in separate page
        # families so an antibody with multiple evidence classes is represented
        # in each view without mixing the claims.
        results = []
        functional_results = []
        negative_results = []
        for antibody_id, interactions in by_target[target_id].items():
            direct = [
                item
                for item in interactions
                if item["relationship"] in DIRECT_POSITIVE_RELATIONSHIPS
            ]
            functional = [
                item
                for item in interactions
                if item["relationship"] in FUNCTIONAL_POSITIVE_RELATIONSHIPS
            ]
            negative = [
                item for item in interactions if _is_negative_relationship(item["relationship"])
            ]
            if direct:
                results.append(make_result(antibody_id, direct))
            if functional:
                functional_results.append(make_result(antibody_id, functional))
            if negative:
                negative_results.append(make_result(antibody_id, negative))

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
            return evidence_score * 100 + len(result["sources"]) * 10

        results.sort(
            key=lambda result: (-result_rank(result), result["antibody"].get("name", "").casefold())
        )
        negative_results.sort(
            key=lambda result: (-result_rank(result), result["antibody"].get("name", "").casefold())
        )
        functional_results.sort(
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

        negative_pages = []
        for page_number, start in enumerate(range(0, len(negative_results), TARGET_PAGE_SIZE), 1):
            page_name = f"negative-page-{page_number:03d}.json"
            page_path = target_dir / page_name
            payload = negative_results[start : start + TARGET_PAGE_SIZE]
            page_path.write_text(
                json.dumps(payload, separators=(",", ":"), ensure_ascii=False), encoding="utf-8"
            )
            page_bytes = page_path.stat().st_size
            target_file_bytes += page_bytes
            largest_target_page_bytes = max(largest_target_page_bytes, page_bytes)
            negative_pages.append({"file": page_name, "count": len(payload), "bytes": page_bytes})

        functional_pages = []
        for page_number, start in enumerate(range(0, len(functional_results), TARGET_PAGE_SIZE), 1):
            page_name = f"functional-page-{page_number:03d}.json"
            page_path = target_dir / page_name
            payload = functional_results[start : start + TARGET_PAGE_SIZE]
            page_path.write_text(
                json.dumps(payload, separators=(",", ":"), ensure_ascii=False), encoding="utf-8"
            )
            page_bytes = page_path.stat().st_size
            target_file_bytes += page_bytes
            largest_target_page_bytes = max(largest_target_page_bytes, page_bytes)
            functional_pages.append({"file": page_name, "count": len(payload), "bytes": page_bytes})

        largest_target_pages = max(
            largest_target_pages, len(pages), len(functional_pages), len(negative_pages)
        )

        def collection_stats(collection: list[dict]) -> dict:
            return {
                "unique_results": len(collection),
                "paired": sum(
                    bool(item["antibody"]["has_heavy"] and item["antibody"]["has_light"])
                    for item in collection
                ),
                "therapeutic": sum(
                    bool(item["antibody"]["therapeutic_status"]) for item in collection
                ),
                "structure_exact": sum(
                    bool(item["antibody"]["exact_structures"]) for item in collection
                ),
                "structure_homologous": sum(
                    bool(item["antibody"]["homologous_structures"]) for item in collection
                ),
            }

        target_stats = {
            **collection_stats(results),
            "functional": len(functional_results),
            "negative": len(negative_results),
            "functional_stats": collection_stats(functional_results),
            "negative_stats": collection_stats(negative_results),
        }
        target_public = {
            "id": target_id,
            "name": target["name"],
            "aliases": sorted(set(target["aliases"]), key=str.casefold),
            "count": target["count"],
            "result_count": len(results),
            "functional_count": len(functional_results),
            "negative_count": len(negative_results),
            "sources": sorted(target["sources"]),
            "relationships": dict(sorted(target["relationships"].items())),
            "dir": target_dir_name,
            "page_count": len(pages),
            "functional_page_count": len(functional_pages),
            "negative_page_count": len(negative_pages),
            "page_size": TARGET_PAGE_SIZE,
            "stats": target_stats,
            "entity": target.get("entity", {}),
            "external_ids": sorted(target.get("external_ids", set())),
            "hierarchy": target.get("hierarchy", {}),
            "children": sorted(target.get("children", set())),
        }
        target_index.append(target_public)
        (target_dir / "index.json").write_text(
            json.dumps(
                {
                    "target": target_public,
                    "pages": pages,
                    "functional_pages": functional_pages,
                    "negative_pages": negative_pages,
                },
                separators=(",", ":"),
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )

    target_index.sort(key=lambda item: (item["name"].casefold(), item["id"]))
    (out / "targets.json").write_text(
        json.dumps(target_index, separators=(",", ":"), ensure_ascii=False), encoding="utf-8"
    )

    antibody_search: dict[str, dict[str, dict]] = defaultdict(dict)
    antibody_exact = build_exact_antibody_index(antibodies)
    antibody_shards: dict[str, dict[str, dict]] = defaultdict(dict)
    sequence_search: dict[str, dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))
    cdr_search = {
        "cdrh3": defaultdict(lambda: defaultdict(set)),
        "cdrl3": defaultdict(lambda: defaultdict(set)),
    }

    for antibody_id, antibody in antibodies.items():
        shard = _antibody_shard(antibody_id)
        search_item = {
            "id": antibody_id,
            "name": antibody.get("name", ""),
            "aliases": antibody.get("aliases", [])[:12],
            "sources": antibody.get("sources", []),
            "direct_targets": antibody.get("direct_targets", [])[:12],
            "functional_targets": antibody.get("functional_targets", [])[:12],
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
            if field in cdr_search:
                cdr_search[field][str(len(sequence_value))][sequence_value].add(antibody_id)

    antibody_search_bytes = 0
    for bucket, payload in antibody_search.items():
        records = sorted(payload.values(), key=lambda item: (item["name"].casefold(), item["id"]))
        path = out / "antibody-search" / f"{bucket}.json"
        path.write_text(
            json.dumps(records, separators=(",", ":"), ensure_ascii=False), encoding="utf-8"
        )
        antibody_search_bytes += path.stat().st_size

    antibody_exact_bytes = 0
    for bucket, payload in antibody_exact.items():
        path = out / "antibody-exact" / f"{bucket}.json"
        path.write_text(
            json.dumps(payload, separators=(",", ":"), ensure_ascii=False), encoding="utf-8"
        )
        antibody_exact_bytes += path.stat().st_size

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

    cdr_stats = {}
    for field, by_length in cdr_search.items():
        directory = out / "sequence" / field
        directory.mkdir(parents=True, exist_ok=True)
        unique = largest = 0
        for length, sequences in by_length.items():
            records = [
                {"sequence": sequence, "antibody_ids": sorted(ids)}
                for sequence, ids in sorted(sequences.items())
            ]
            path = directory / f"{int(length):02d}.json"
            path.write_text(
                json.dumps(records, separators=(",", ":"), ensure_ascii=False), encoding="utf-8"
            )
            unique += len(records)
            largest = max(largest, path.stat().st_size)
        cdr_stats[field] = {
            "unique_sequences": unique,
            "bucket_files": len(by_length),
            "largest_bucket_bytes": largest,
        }

    similarity_stats = write_similarity_indexes(antibodies, out)
    cluster_stats = write_cluster_indexes(antibodies, out)

    review = _review_candidates(targets, antibody_direct_targets)
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

    if validate_examples:
        examples = json.loads((CONFIG / "search_examples.json").read_text(encoding="utf-8"))
        target_lookup = {}
        for target in target_index:
            for value in [
                target.get("id", ""),
                target.get("name", ""),
                *(target.get("aliases", []) or []),
                *(target.get("external_ids", []) or []),
            ]:
                normalized = normalize_search_term(value)
                if normalized:
                    target_lookup.setdefault(normalized, set()).add(target["id"])
        for query in examples.get("target_examples", []):
            normalized = normalize_search_term(query)
            if not target_lookup.get(normalized):
                raise ValueError(f"advertised target example has no exact generated match: {query!r}")
        for query in examples.get("antibody_examples", []):
            normalized = normalize_search_term(query)
            bucket = antibody_exact.get(search_bucket(normalized), {})
            if not bucket.get(normalized):
                raise ValueError(
                    f"advertised antibody example has no exact generated match: {query!r}"
                )

    exact_terms = sum(len(payload) for payload in antibody_exact.values())
    exact_aliases = sum(
        len(
            {
                normalize_search_term(term)
                for term in antibody.get("aliases", [])
                if normalize_search_term(term)
            }
        )
        for antibody in antibodies.values()
    )
    multi_entity_terms = sum(
        1
        for payload in antibody_exact.values()
        for ids in payload.values()
        if len(ids) > 1
    )
    return {
        "antibodies": len(antibodies),
        "targets": len(targets),
        "interactions": interaction_count,
        "indexed_interactions": indexed_interaction_count,
        "literature_mentions": literature_mention_count,
        "source_counts": dict(source_counts),
        "files": {
            "targets_index_bytes": (out / "targets.json").stat().st_size,
            "target_page_bytes": target_file_bytes,
            "largest_target_page_bytes": largest_target_page_bytes,
            "largest_target_page_count": largest_target_pages,
            "antibody_search_index_bytes": antibody_search_bytes,
            "antibody_exact_index_bytes": antibody_exact_bytes,
            "sequence_search_index_bytes": sequence_search_bytes,
        },
        "antibody_exact_terms": exact_terms,
        "antibody_exact_aliases": exact_aliases,
        "antibody_exact_multi_entity_terms": multi_entity_terms,
        "target_review_candidates": len(review),
        "cdr_indexes": cdr_stats,
        "similarity_indexes": similarity_stats,
        "cluster_indexes": cluster_stats,
    }


def compile_data(
    source_paths: dict[str, Path],
    source_config: dict,
    out: Path,
    max_records: int | None = None,
    validate_examples: bool = False,
) -> dict:
    antibodies: dict[str, dict] = {}
    raw_interactions: list[dict] = []
    source_counts: dict[str, dict] = defaultdict(lambda: {"records": 0, "interactions": 0})
    quarantined_records = 0
    multispecific_construct_ids: set[str] = set()
    unassigned_multispecific_arm_ids: set[str] = set()

    for source, path in source_paths.items():
        config = source_config[source]
        adapter = ADAPTERS[config["adapter"]]
        source_record_groups: set[str] = set()
        for index, (observation, interactions) in enumerate(
            adapter(path, config.get("homepage", ""))
        ):
            record_group = observation.construct.get("id") or (
                f"{source}:{observation.record_id or index}"
            )
            is_new_record = record_group not in source_record_groups
            if (
                max_records is not None
                and is_new_record
                and len(source_record_groups) >= max_records
            ):
                break
            source_record_groups.add(record_group)
            antibody_id = observation.identity()
            antibody = antibodies.setdefault(
                antibody_id, {"id": antibody_id, "name": observation.name}
            )
            merge_antibody(antibody, observation)
            if is_new_record:
                source_counts[source]["records"] += 1
            if observation.metadata.get("multispecific") and observation.construct.get("id"):
                multispecific_construct_ids.add(observation.construct["id"])
            if (
                observation.arm.get("id")
                and observation.arm.get("target_assignment_status") == "unavailable_no_arm_mapping"
            ):
                unassigned_multispecific_arm_ids.add(observation.arm["id"])
            if observation.metadata.get("sequence_quarantine"):
                quarantined_records += 1
                # Preserve the observation in the antibody/provenance index,
                # but quarantine all target claims until its sequence is
                # corrected upstream.
                continue
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
                        "assertion_origin": interaction.assertion_origin,
                        "source_url": observation.source_url,
                        "record_url": interaction.record_url or observation.record_url,
                        "link_scope": (
                            interaction.link_scope
                            if interaction.record_url
                            else observation.link_scope
                        ),
                        "measurements": interaction.measurements,
                        "target_external_id": interaction.target_external_id,
                        "assay_ids": interaction.assay_ids,
                        "receptor_group_id": interaction.receptor_group_id,
                    }
                )

    stats = write_indexes(
        antibodies,
        raw_interactions,
        out,
        source_counts,
        validate_examples=validate_examples,
    )
    stats["quarantined_records"] = quarantined_records
    stats["multispecific_constructs"] = len(multispecific_construct_ids)
    stats["unassigned_multispecific_arms"] = len(unassigned_multispecific_arm_ids)
    return stats


def rebuild_similarity_only(out: Path) -> dict:
    antibodies = {}
    for shard_path in sorted((out / "antibodies").glob("*.json")):
        antibodies.update(json.loads(shard_path.read_text(encoding="utf-8")))
    if not antibodies:
        raise FileNotFoundError(f"no antibody shards found under {out}")
    stats = write_similarity_indexes(antibodies, out)
    clusters = write_cluster_indexes(antibodies, out)
    manifest_path = out / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.setdefault("stats", {})["similarity_indexes"] = stats
    manifest["stats"]["cluster_indexes"] = clusters
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return {"similarity_indexes": stats, "cluster_indexes": clusters}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Build the static PAIRS dataset")
    parser.add_argument("--output", default=str(ROOT / "data" / DATA_SUBDIR))
    parser.add_argument("--cache", default=str(ROOT / ".cache" / "sources"))
    parser.add_argument(
        "--offline", action="store_true", help="Use source files already present in --cache"
    )
    parser.add_argument(
        "--similarity-only",
        action="store_true",
        help="Rebuild similarity and sequence-cluster indexes from existing antibody shards",
    )
    parser.add_argument(
        "--max-records",
        type=int,
        default=None,
        help="Limit records per source (for local/demo builds)",
    )
    parser.add_argument(
        "--source",
        action="append",
        dest="sources",
        help="Build only selected source key; repeatable",
    )
    parser.add_argument(
        "--allow-partial", action="store_true", help="Continue when a source download fails"
    )
    args = parser.parse_args(argv)

    if args.similarity_only:
        print(json.dumps(rebuild_similarity_only(Path(args.output)), indent=2))
        return 0

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
                if source_config.get("pagination") == "postgrest_csv":
                    info = download_postgrest_csv(resolved_url, path)
                else:
                    info = download(
                        resolved_url,
                        path,
                        accept=source_config.get("accept", ""),
                    )
                info.update(discovery)
                info["download_url"] = resolved_url
                if source == "iedb":
                    info["support"] = download_iedb_support(cache)
            elif not path.exists():
                raise FileNotFoundError(path)
            elif source == "iedb":
                for support_name in ("iedb-bcr-search.csv", "iedb-bcell-export.csv"):
                    if not (cache / support_name).exists():
                        raise FileNotFoundError(cache / support_name)
                info["support"] = inspect_cached_iedb_support(cache)
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
    stats = compile_data(
        source_paths,
        config,
        out,
        args.max_records,
        validate_examples=args.max_records is None and not args.allow_partial,
    )
    for source, counts in stats.get("source_counts", {}).items():
        statuses.setdefault(source, {}).update(counts)

    snapshot = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "data_contract_revision": DATA_CONTRACT_REVISION,
        "app_version": APP_VERSION,
        "data_path": f"data/{DATA_SUBDIR}",
        "snapshot": snapshot,
        "snapshot_date": snapshot,
        "antibody_count": stats["antibodies"],
        "target_count": stats["targets"],
        "interaction_count": stats["interactions"],
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
