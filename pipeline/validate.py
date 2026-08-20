from __future__ import annotations

import json
import math
import sys
from pathlib import Path

from .build import DATA_CONTRACT_REVISION, _sequence_sha, _sequence_signature, sequence_quality
from .model import SequenceNormalizationError, sequence as normalize_sequence, sequence_contract
from .search import (
    antibody_search_terms,
    iter_exact_antibody_terms,
    normalize_search_term,
    search_bucket,
)

EXPECTED_SCHEMA = 4


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
    if manifest.get("data_contract_revision") != DATA_CONTRACT_REVISION:
        errors.append("unexpected data_contract_revision")
    if manifest.get("stats", {}).get("targets") != len(targets):
        errors.append("target count mismatch")
    if not antibody_shards:
        errors.append("no antibody shards")

    antibody_ids: set[str] = set()
    antibodies_by_id: dict[str, dict] = {}
    unassigned_construct_records: set[tuple[str, str]] = set()
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
            antibody = payload[antibody_id]
            if "direct_targets" not in antibody or not isinstance(
                antibody.get("direct_targets"), list
            ):
                errors.append(f"V4 antibody lacks direct_targets list: {antibody_id}")
            if antibody.get("sequence_quality") and antibody[
                "sequence_quality"
            ] != sequence_quality(antibody):
                errors.append(f"inconsistent sequence quality metadata on {antibody_id}")
            for field in ["vh_nt_source", "vl_nt_source"]:
                source_nt = antibody.get(field, "")
                if source_nt and (set(source_nt) - set("ACGTRYSWKMBDHVN")):
                    errors.append(f"invalid source nucleotide sequence on {antibody_id}: {field}")
                if source_nt and not antibody.get("nucleotide_provenance", {}).get(
                    "VH" if field.startswith("vh_") else "VL"
                ):
                    errors.append(
                        f"source nucleotide sequence lacks provenance on {antibody_id}: {field}"
                    )
                if source_nt:
                    provenance = antibody.get("nucleotide_provenance", {}).get(
                        "VH" if field.startswith("vh_") else "VL", {}
                    )
                    if provenance.get("scope") not in {
                        "variable_domain",
                        "full_length_bcr",
                        "unknown",
                    }:
                        errors.append(
                            f"source nucleotide sequence has invalid scope on {antibody_id}: {field}"
                        )
            source_records = {
                (record.get("source", ""), record.get("record_id", ""))
                for record in antibody.get("source_records", [])
            }
            for nucleotide in antibody.get("source_nucleotide_records", []):
                key = (nucleotide.get("source", ""), nucleotide.get("source_record_id", ""))
                if key not in source_records:
                    errors.append(
                        f"source nucleotide lacks source record on {antibody_id}: {key}"
                    )
                if nucleotide.get("chain") not in {"VH", "VL"}:
                    errors.append(f"invalid source nucleotide chain on {antibody_id}")
                if nucleotide.get("scope") not in {
                    "variable_domain",
                    "full_length_bcr",
                    "unknown",
                }:
                    errors.append(f"invalid source nucleotide scope on {antibody_id}")
                if nucleotide.get("sequence") and set(nucleotide["sequence"]) - set(
                    "ACGTRYSWKMBDHVN"
                ):
                    errors.append(f"invalid source nucleotide record on {antibody_id}")
            for source_record in antibody.get("source_records", []):
                record_key = (source_record.get("source", ""), source_record.get("record_id", ""))
                for nucleotide in source_record.get("nucleotide_records", []):
                    nucleotide_key = (
                        nucleotide.get("source", "") or record_key[0],
                        nucleotide.get("source_record_id", "") or record_key[1],
                    )
                    if nucleotide_key != record_key:
                        errors.append(
                            f"source record nucleotide provenance mismatch on {antibody_id}"
                        )
            annotations = antibody.get("chain_annotations", {})
            if not isinstance(annotations, dict) or any(
                chain not in {"VH", "VL"} or not isinstance(annotation, dict)
                for chain, annotation in annotations.items()
            ):
                errors.append(f"invalid chain annotations on {antibody_id}")
            for source_record in antibody.get("source_records", []):
                if source_record.get("link_scope") == "record" and not source_record.get(
                    "record_url"
                ):
                    errors.append(f"record-scoped provenance has no URL on {antibody_id}")
            constructs = antibody.get("constructs", [])
            if antibody.get("metadata", {}).get("multispecific") and (
                not constructs or not antibody.get("arms")
            ):
                errors.append(f"multispecific antibody lacks construct/arm context: {antibody_id}")
            for construct in constructs:
                if construct.get("target_assignment_status") == "unavailable_no_arm_mapping":
                    unassigned_construct_records.add(
                        (construct.get("source", ""), construct.get("source_record_id", ""))
                    )

    if manifest.get("stats", {}).get("antibodies") != len(antibody_ids):
        errors.append("antibody count mismatch")

    target_ids: set[str] = set()
    referenced_antibodies: set[str] = set()
    counted_results = 0
    measurement_keys: set[tuple] = set()
    for target in targets:
        target_id = target.get("id", "")
        if target_id in target_ids:
            errors.append(f"duplicate target id {target_id}")
        target_ids.add(target_id)
        entity = target.get("entity", {})
        if entity and not entity.get("mapping_provenance"):
            errors.append(f"target entity lacks mapping provenance: {target_id}")
        target_dir = data_dir / "targets" / target.get("dir", "")
        index_path = target_dir / "index.json"
        if not index_path.exists():
            errors.append(f"missing target index for {target_id}")
            continue
        target_manifest = _read_json(index_path)
        target_result_count = 0

        def check_pages(
            manifest_key: str,
            count_key: str,
            allowed_relationships: tuple[str, ...],
            *,
            prefixes: bool = False,
        ) -> int:
            pages = target_manifest.get(manifest_key, [])
            if len(pages) != target.get(count_key, 0):
                errors.append(f"{manifest_key} count mismatch for {target_id}")
            row_count = 0
            for page in pages:
                page_path = target_dir / page.get("file", "")
                if not page_path.exists():
                    errors.append(f"missing target page {target_id}/{page.get('file')}")
                    continue
                rows = _read_json(page_path)
                if len(rows) != page.get("count"):
                    errors.append(f"page row count mismatch for {target_id}/{page.get('file')}")
                row_count += len(rows)
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
                        relationship = interaction.get("relationship", "")
                        permitted = (
                            relationship.startswith(allowed_relationships)
                            if prefixes
                            else relationship in allowed_relationships
                        )
                        if not permitted:
                            errors.append(
                                f"unexpected relationship in {manifest_key} for {target_id}: "
                                f"{interaction.get('relationship')}"
                            )
                        if (
                            interaction.get("source", ""),
                            interaction.get("source_record_id", ""),
                        ) in unassigned_construct_records:
                            errors.append(
                                "unassigned multispecific construct created an arm target claim: "
                                f"{interaction.get('source_record_id')}"
                            )
                        for measurement in interaction.get("measurements", []):
                            if measurement.get("metric") not in {
                                "KD",
                                "KON",
                                "KOFF",
                                "IC50",
                                "EC50",
                            }:
                                errors.append(
                                    f"invalid measurement metric on {interaction.get('id')}"
                                )
                            if measurement.get("qualifier", "") not in {
                                "",
                                "<",
                                "<=",
                                "=",
                                ">=",
                                ">",
                                "~",
                            }:
                                errors.append(
                                    f"invalid measurement qualifier on {interaction.get('id')}"
                                )
                            value = measurement.get("value")
                            if value is None and not measurement.get("raw_value"):
                                errors.append(f"empty measurement on {interaction.get('id')}")
                            if value is not None and not isinstance(value, (int, float)):
                                errors.append(
                                    f"nonnumeric measurement value on {interaction.get('id')}"
                                )
                            if isinstance(value, (int, float)) and not math.isfinite(value):
                                errors.append(
                                    f"nonfinite measurement value on {interaction.get('id')}"
                                )
                            if not interaction.get("source_record_id"):
                                errors.append(f"orphan measurement on {interaction.get('id')}")
                            assay_id = str(measurement.get("assay_id", ""))
                            if assay_id and assay_id not in {
                                str(item) for item in interaction.get("assay_ids", [])
                            }:
                                errors.append(
                                    f"measurement assay linkage mismatch on {interaction.get('id')}"
                                )
                            measurement_key = (
                                antibody_id,
                                interaction.get("source"),
                                interaction.get("receptor_group_id"),
                                assay_id,
                                measurement.get("metric"),
                                measurement.get("raw_value"),
                                measurement.get("unit"),
                                measurement.get("qualifier"),
                                interaction.get("target_external_id"),
                            )
                            if measurement_key in measurement_keys:
                                errors.append(
                                    f"duplicate measurement linkage on {interaction.get('id')}"
                                )
                            measurement_keys.add(measurement_key)
            return row_count

        target_result_count = check_pages("pages", "page_count", ("binds", "targets"))
        functional_result_count = check_pages(
            "functional_pages", "functional_page_count", ("neutralizes", "protects")
        )
        negative_result_count = check_pages(
            "negative_pages",
            "negative_page_count",
            ("does_not_", "not_"),
            prefixes=True,
        )
        if target_result_count != target.get("result_count"):
            errors.append(f"target result count mismatch for {target_id}")
        if functional_result_count != target.get("functional_count", 0):
            errors.append(f"target functional count mismatch for {target_id}")
        if negative_result_count != target.get("negative_count", 0):
            errors.append(f"target negative count mismatch for {target_id}")
        counted_results += target_result_count

    hierarchy = {
        target.get("id", ""): target.get("hierarchy", {}).get("parent_id", "")
        for target in targets
        if target.get("hierarchy", {}).get("parent_id")
    }
    for child_id, parent_id in hierarchy.items():
        if parent_id not in target_ids:
            errors.append(f"target hierarchy references missing parent: {child_id} -> {parent_id}")
        seen = {child_id}
        cursor = parent_id
        while cursor:
            if cursor in seen:
                errors.append(f"target hierarchy cycle at {child_id}")
                break
            seen.add(cursor)
            cursor = hierarchy.get(cursor, "")
    for target in targets:
        for child_id in target.get("children", []):
            if hierarchy.get(child_id) != target.get("id"):
                errors.append(
                    f"target child edge is not reciprocal: {target.get('id')} -> {child_id}"
                )

    if not any((data_dir / "antibody-search").glob("*.json")):
        errors.append("no antibody search shards")
    exact_directory = data_dir / "antibody-exact"
    exact_lookup: dict[str, list[str]] = {}
    if not exact_directory.exists() or not any(exact_directory.glob("*.json")):
        errors.append("no antibody exact-name shards")
    else:
        for bucket_path in exact_directory.glob("*.json"):
            payload = _read_json(bucket_path)
            if not isinstance(payload, dict):
                errors.append(f"antibody exact shard is not an object: {bucket_path.name}")
                continue
            for term, ids in payload.items():
                if not isinstance(term, str) or not term or normalize_search_term(term) != term:
                    errors.append(
                        f"invalid normalized antibody exact term: {bucket_path.name}/{term!r}"
                    )
                    continue
                if search_bucket(term) != bucket_path.stem:
                    errors.append(f"antibody exact term in wrong shard: {term}")
                if not isinstance(ids, list) or not ids:
                    errors.append(f"antibody exact term has no IDs: {term}")
                    continue
                if any(
                    not isinstance(antibody_id, str) or not antibody_id for antibody_id in ids
                ):
                    errors.append(f"antibody exact term has malformed ID: {term}")
                if len(ids) != len(set(ids)):
                    errors.append(f"duplicate IDs in antibody exact posting: {term}")
                exact_lookup[term] = ids
                for antibody_id in ids:
                    if antibody_id not in antibody_ids:
                        errors.append(
                            f"antibody exact index {bucket_path.name} references missing antibody {antibody_id}"
                        )
                    antibody = antibodies_by_id.get(antibody_id)
                    if antibody and not any(
                        normalize_search_term(value) == term
                        for value in antibody_search_terms(antibody)
                    ):
                        errors.append(
                            f"antibody exact term is not present on referenced antibody {term}: {antibody_id}"
                        )

        for antibody_id, raw_term, normalized in iter_exact_antibody_terms(antibodies_by_id):
            if antibody_id not in exact_lookup.get(normalized, []):
                errors.append(
                    f"antibody exact round-trip failed for {antibody_id}: {raw_term!r} -> {normalized!r}"
                )
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
                    continue
                field = match.get("field", "")
                if field not in {"heavy", "light", "cdrh3", "cdrl3"}:
                    errors.append(
                        f"sequence index {sequence_path.name} has unsupported field {field!r}"
                    )
                    continue
                indexed = antibodies_by_id[antibody_id].get(field, "")
                try:
                    if normalize_sequence(indexed) != indexed:
                        errors.append(
                            f"sequence index value is not normalized on {antibody_id}: {field}"
                        )
                except SequenceNormalizationError:
                    errors.append(f"sequence index value is invalid on {antibody_id}: {field}")
                if _sequence_sha(indexed) != sequence_hash:
                    errors.append(
                        f"sequence index hash does not match {antibody_id}: {field}"
                    )
                if match.get("length") != len(indexed):
                    errors.append(f"sequence index length does not match {antibody_id}: {field}")

    contract = sequence_contract()
    if contract.get("version") != 1 or not contract.get("exact_alphabet"):
        errors.append("unsupported sequence normalization contract")

    for field in ("cdrh3", "cdrl3"):
        directory = data_dir / "sequence" / field
        if not directory.exists():
            errors.append(f"missing {field} index directory")
            continue
        for bucket_path in directory.glob("*.json"):
            try:
                bucket_length = int(bucket_path.stem)
            except ValueError:
                errors.append(f"invalid {field} bucket name: {bucket_path.name}")
                continue
            for record in _read_json(bucket_path):
                sequence = record.get("sequence", "")
                if len(sequence) != bucket_length:
                    errors.append(f"{field} sequence in wrong bucket: {sequence}")
                for antibody_id in record.get("antibody_ids", []):
                    if antibody_id not in antibody_ids:
                        errors.append(f"{field} index references missing antibody {antibody_id}")
                    elif sequence != antibodies_by_id[antibody_id].get(field, ""):
                        errors.append(
                            f"{field} index sequence not present on antibody {antibody_id}"
                        )

    for field in ("heavy", "light"):
        directory = data_dir / "similarity" / field
        metadata_path = directory / "index.json"
        if not metadata_path.exists():
            errors.append(f"missing {field} similarity index metadata")
            continue
        metadata = _read_json(metadata_path)
        if metadata.get("version") != 1 or metadata.get("k") != 5:
            errors.append(f"unsupported {field} similarity index contract")
        signatures = {
            antibody_id: set(_sequence_signature(antibody.get(field, "")))
            for antibody_id, antibody in antibodies_by_id.items()
        }
        for bucket_path in directory.glob("[0-9a-f].json"):
            for signature_hash, ids in _read_json(bucket_path).items():
                if not signature_hash.startswith(bucket_path.stem):
                    errors.append(f"similarity hash in wrong bucket: {signature_hash}")
                if len(ids) != len(set(ids)):
                    errors.append(f"duplicate IDs in similarity posting: {signature_hash}")
                for antibody_id in ids:
                    if antibody_id not in antibody_ids:
                        errors.append(f"similarity index references missing antibody {antibody_id}")
                    elif signature_hash not in signatures[antibody_id]:
                        errors.append(
                            f"{field} similarity signature absent from antibody {antibody_id}"
                        )

    cluster_metadata_path = data_dir / "clusters" / "index.json"
    if not cluster_metadata_path.exists():
        errors.append("missing sequence-cluster metadata")
    else:
        cluster_metadata = _read_json(cluster_metadata_path)
        if (
            cluster_metadata.get("version") != 1
            or cluster_metadata.get("minimum_coverage") != 90
            or cluster_metadata.get("thresholds") != [99, 95, 90]
        ):
            errors.append("unsupported sequence-cluster contract")
        for scope in ("heavy", "light", "paired"):
            required_fields = ("heavy", "light") if scope == "paired" else (scope,)
            for threshold in (99, 95, 90):
                directory = data_dir / "clusters" / scope / str(threshold)
                index_path = directory / "index.json"
                if not index_path.exists():
                    errors.append(f"missing {scope} {threshold}% cluster index")
                    continue
                index = _read_json(index_path)
                lookup = {}
                for lookup_path in (directory / "lookup").glob("*.json"):
                    lookup.update(_read_json(lookup_path))
                seen_members = set()
                indexed_lookup = {}
                for cluster in index.get("clusters", []):
                    members = cluster.get("members", [])
                    if len(members) < 2 or members != sorted(set(members)):
                        errors.append(f"invalid membership in cluster {cluster.get('id')}")
                    if cluster.get("representative_id") != min(members, default=""):
                        errors.append(f"nondeterministic representative in {cluster.get('id')}")
                    for antibody_id in members:
                        if antibody_id in seen_members:
                            errors.append(
                                f"overlapping {scope} {threshold}% clusters: {antibody_id}"
                            )
                        seen_members.add(antibody_id)
                        antibody = antibodies_by_id.get(antibody_id)
                        if not antibody or not all(
                            antibody.get(field) for field in required_fields
                        ):
                            errors.append(
                                f"{scope} cluster references incompatible antibody {antibody_id}"
                            )
                        indexed_lookup[antibody_id] = cluster.get("id")
                if lookup != indexed_lookup:
                    errors.append(f"{scope} {threshold}% cluster lookup mismatch")

    stats = manifest.get("stats", {})
    iedb_source = manifest.get("sources", {}).get("iedb", {})
    if iedb_source.get("ok"):
        support = iedb_source.get("support", {})
        if (
            support.get("linked_assay_ids", 0) <= 0
            or support.get("bcell_export_rows", 0) <= 0
            or support.get("bcr_search", {}).get("rows", 0) <= 0
        ):
            errors.append("IEDB support linkage is empty")
        if not measurement_keys:
            errors.append("IEDB is included but no linked measurements were emitted")
    if stats.get("interactions", 0) <= 0:
        errors.append("interaction count is zero")
    if targets and not referenced_antibodies:
        errors.append("target pages reference no antibodies")
    if counted_results <= 0:
        errors.append("target result pages are empty")
    return errors


def main() -> int:
    data_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "data/v4")
    errors = validate(data_dir)
    if errors:
        print("\n".join("ERROR: " + error for error in errors), file=sys.stderr)
        return 1
    print(f"Validation OK: {data_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
