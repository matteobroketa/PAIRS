from __future__ import annotations

import csv
import gzip
import io
import json
import re
from collections import defaultdict
from urllib.parse import quote
from pathlib import Path

from .model import (
    AntibodyObservation,
    InteractionObservation,
    NA,
    SequenceNormalizationError,
    digest,
    infer_sequence_architecture,
    measurement,
    nucleotide_sequence,
    sequence,
    split_values,
    text,
)

csv.field_size_limit(50_000_000)

THERA_SUMMARY_URL = (
    "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/therasabdab/therasummary/?INN="
)


def _open_csv(path: Path):
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    decoded = raw.decode("utf-8-sig", errors="replace")
    return csv.DictReader(io.StringIO(decoded, newline=""))


def _open_json(path: Path):
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return json.loads(raw.decode("utf-8-sig"))


def _require(reader, fields: list[str], source: str):
    available = set(reader.fieldnames or [])
    missing = [f for f in fields if f not in available]
    if missing:
        raise ValueError(f"{source} schema changed; missing columns: {', '.join(missing)}")
    return reader


def _source_accessions(row: dict) -> list[dict]:
    """Extract only antibody-chain accessions explicitly supplied by a source."""
    output = []
    for field, value in row.items():
        field_text = text(field)
        field_key = field_text.casefold().replace("_", " ")
        if not value or not (
            "accession" in field_key
            or "genbank" in field_key
            or "ena" in field_key
            or "embl" in field_key
        ):
            continue
        if not re.search(r"(?:vh|heavy|variable heavy|vl|light|variable light)", field_key):
            continue
        chain = "VH" if re.search(r"(?:vh|heavy|variable heavy)", field_key) else "VL"
        database = (
            "ena"
            if "ena" in field_key
            else "embl"
            if "embl" in field_key
            else "genbank"
        )
        scope = "variable_domain" if "variable" in field_key else "unknown"
        for raw_accession in split_values(value, separators=r"[;,|\s]+"):
            accession = raw_accession.strip()
            if not re.fullmatch(r"[A-Za-z]{1,5}[A-Za-z0-9_]*\d+(?:\.\d+)?", accession):
                continue
            output.append(
                {
                    "accession": accession,
                    "database": database,
                    "chain": chain,
                    "scope": scope,
                    "source_field": field_text,
                }
            )
    unique = {}
    for item in output:
        unique[(item["database"], item["accession"], item["chain"])] = item
    return list(unique.values())


def _structures(value: str) -> list[str]:
    s = text(value)
    if not s:
        return []
    hits = re.findall(r"\b[0-9][A-Za-z0-9]{3}\b", s)
    return sorted(set(x.upper() for x in hits))


def _thera_record_url(therapeutic: str) -> str:
    """Build the documented Thera-SAbDab summary URL for one INN."""
    return THERA_SUMMARY_URL + quote(therapeutic, safe="")


def _thera_structure_tiers(row: dict) -> dict[str, list[str]]:
    """Keep Thera-SAbDab's exact and homologous structure matches distinct."""
    return {
        "100%": _structures(row.get("100% SI Structure")),
        "99%": _structures(row.get("99% SI Structure")),
        "95-98%": _structures(row.get("95-98% SI Structure")),
    }


def _thera_multispecific(row: dict) -> bool:
    """Return whether a Thera-SAbDab row cannot be mapped arm-by-arm yet.

    Thera-SAbDab keeps the construct's targets in one column, while the
    additional ``*Sequence(ifbispec)`` columns contain a second variable
    domain set.  Until those domains are represented as separate arms, any
    target emitted for the row would be an unsafe sequence-to-target claim.
    Treat an explicit multi-specific format, populated auxiliary arm fields,
    or a semicolon-delimited target list as ambiguous and quarantine it.
    """
    format_value = text(row.get("Format"))
    if re.search(r"\b(?:bi|tri|tetra|multi|quadri|penta)spec", format_value, re.IGNORECASE):
        return True

    if any(
        text(row.get(field)) for field in ("HeavySequence(ifbispec)", "LightSequence(ifbispec)")
    ):
        return True

    # A semicolon is the source convention for separate bispecific target
    # groups.  Keep ambiguous rows out of target retrieval even when Format
    # is missing or stale in the download.
    return len(split_values(row.get("Target"))) > 1


def _thera_arm_sequence(raw_heavy: object, raw_light: object) -> tuple[str, str, bool]:
    """Normalize one candidate Thera arm and report whether it is usable."""
    try:
        heavy = sequence(raw_heavy)
        light = sequence(raw_light)
    except SequenceNormalizationError:
        return "", "", False
    return heavy, light, bool(heavy or light)


def _thera_arm_mapping(row: dict, target_groups: list[str]) -> tuple[list[dict], str]:
    """Build arm records, mapping targets only for an unambiguous two-arm row.

    The downloadable schema has one primary VH/VL pair and optional
    ``(ifbispec)`` VH/VL fields.  It does *not* provide an arm-specific target
    column, so the target list cannot safely be assigned by position.  The
    explicit-column branch is intentionally future-facing: only a source
    field naming the target for a particular arm may enable interactions.
    """
    candidates = [
        ("primary", row.get("HeavySequence"), row.get("LightSequence")),
        ("ifbispec", row.get("HeavySequence(ifbispec)"), row.get("LightSequence(ifbispec)")),
    ]
    arms: list[dict] = []
    for role, raw_heavy, raw_light in candidates:
        heavy, light, available = _thera_arm_sequence(raw_heavy, raw_light)
        if available:
            arms.append({"role": role, "heavy": heavy, "light": light})

    explicit_targets = {
        "primary": text(row.get("Target (HeavySequence)")) or text(row.get("Target (primary)")),
        "ifbispec": text(row.get("Target (HeavySequence(ifbispec))"))
        or text(row.get("Target (ifbispec)")),
    }
    arm_mapped = bool(arms) and all(explicit_targets.get(arm["role"]) for arm in arms)
    status = "arm_mapped_explicit_source_field" if arm_mapped else "unavailable_no_arm_mapping"
    for index, arm in enumerate(arms, 1):
        arm["label"] = f"arm_{index}"
        if arm_mapped:
            arm["target_raw"] = explicit_targets[arm["role"]]
        arm["target_assignment_status"] = status
    return arms, status


def _obs(source, rid, name, heavy="", light="", **kw):
    metadata = dict(kw.pop("metadata", {}) or {})
    explicit_architecture = kw.pop("sequence_architecture", "")
    sequence_accessions = list(kw.pop("sequence_accessions", []) or [])
    invalid_sequences = {}
    normalized = {}
    for field, value in {
        "heavy": heavy,
        "light": light,
        "cdrh1": kw.pop("cdrh1", ""),
        "cdrh2": kw.pop("cdrh2", ""),
        "cdrh3": kw.pop("cdrh3", ""),
        "cdrl1": kw.pop("cdrl1", ""),
        "cdrl2": kw.pop("cdrl2", ""),
        "cdrl3": kw.pop("cdrl3", ""),
    }.items():
        raw = text(value)
        try:
            normalized[field] = sequence(raw)
        except SequenceNormalizationError as exc:
            # Keep the record available for provenance/name search, but do
            # not let an invalid source value become a target-bearing sequence.
            normalized[field] = ""
            invalid_sequences[field] = {"raw": raw, "error": str(exc)}
    normalized_nt = {}
    for field in ["vh_nt_source", "vl_nt_source"]:
        raw = text(kw.pop(field, ""))
        try:
            normalized_nt[field] = nucleotide_sequence(raw)
        except SequenceNormalizationError as exc:
            normalized_nt[field] = ""
            invalid_sequences[field] = {"raw": raw, "error": str(exc)}
    nucleotide_provenance = dict(kw.pop("nucleotide_provenance", {}) or {})
    for field, chain in [("vh_nt_source", "VH"), ("vl_nt_source", "VL")]:
        if normalized_nt[field]:
            nucleotide_provenance.setdefault(
                chain,
                {"source": source, "source_record_id": text(rid) or text(name)},
            )
    source_nucleotide_records = []
    for field, chain in [("vh_nt_source", "VH"), ("vl_nt_source", "VL")]:
        if not normalized_nt[field]:
            continue
        provenance = nucleotide_provenance.get(chain, {})
        source_nucleotide_records.append(
            {
                "source": provenance.get("source") or source,
                "source_record_id": provenance.get("source_record_id")
                or text(rid)
                or text(name),
                "chain": chain,
                "sequence": normalized_nt[field],
                "scope": provenance.get("scope") or "unknown",
                "source_field": provenance.get("source_field") or field,
            }
        )
    if invalid_sequences:
        metadata.update(
            {
                "sequence_quarantine": True,
                "sequence_quarantine_fields": invalid_sequences,
            }
        )
    architecture = infer_sequence_architecture(
        normalized["heavy"],
        normalized["light"],
        tuple(
            normalized[field]
            for field in ("cdrh1", "cdrh2", "cdrh3", "cdrl1", "cdrl2", "cdrl3")
        ),
        normalized_nt["vh_nt_source"],
        normalized_nt["vl_nt_source"],
        kw.get("format", ""),
        explicit_architecture,
    )
    metadata.setdefault("sequence_architecture", architecture)
    return AntibodyObservation(
        source=source,
        record_id=text(rid) or text(name),
        name=text(name) or text(rid) or "Unnamed antibody",
        heavy=normalized["heavy"],
        light=normalized["light"],
        cdrh3=normalized["cdrh3"],
        cdrl3=normalized["cdrl3"],
        cdrh1=normalized["cdrh1"],
        cdrh2=normalized["cdrh2"],
        cdrl1=normalized["cdrl1"],
        cdrl2=normalized["cdrl2"],
        vh_nt_source=normalized_nt["vh_nt_source"],
        vl_nt_source=normalized_nt["vl_nt_source"],
        sequence_architecture=architecture,
        sequence_accessions=sequence_accessions,
        nucleotide_provenance=nucleotide_provenance,
        source_nucleotide_records=source_nucleotide_records,
        metadata=metadata,
        **kw,
    )


def plabdab(path: Path, source_url: str = ""):
    reader = _require(
        _open_csv(path), ["ID", "heavy_sequence", "light_sequence", "targets_mentioned"], "PLAbDab"
    )
    for row in reader:
        target_group = text(row.get("targets_mentioned"))
        ab = _obs(
            "plabdab",
            row.get("ID"),
            row.get("ID"),
            row.get("heavy_sequence"),
            row.get("light_sequence"),
            sequence_accessions=_source_accessions(row),
            organism=text(row.get("organism")),
            format="paired antibody" if text(row.get("light_sequence")) else "antibody",
            reference=text(row.get("reference_title")),
            source_url=source_url,
            metadata={
                "pairing": text(row.get("pairing")),
                "update_date": text(row.get("update_date")),
            },
        )
        aid = ab.identity()
        interactions = []
        for mentioned_target in split_values(target_group, separators=r";"):
            interactions.append(
                InteractionObservation(
                    antibody_id=aid,
                    source="plabdab",
                    source_record_id=ab.record_id,
                    target_raw=mentioned_target,
                    relationship="mentioned_with",
                    evidence="LITERATURE_METADATA",
                    reference=ab.reference,
                    note=(
                        "Target term extracted from a PLAbDab literature mention group. "
                        "Co-mentioned terms are indexed independently and are not treated as synonyms "
                        "or as direct binding measurements."
                    ),
                )
            )
        yield ab, interactions


def _sabdab_segment_sequence(segment: dict | None) -> str:
    """Recover an unnumbered variable-domain sequence from SAbDab2 JSON."""
    if not segment:
        return ""
    direct = text(segment.get("sequence"))
    if direct:
        return direct.replace("-", "")
    numbering = segment.get("numbering") or []
    return "".join(
        residue
        for _, residue in numbering
        if text(residue) and text(residue) != "-"
    )


def _sabdab_pdb_id(value: object) -> str:
    raw = text(value)
    match = re.search(r"([0-9A-Za-z]{4})$", raw)
    return match.group(1).upper() if match else raw.upper()


def _sabdab_labels(value: object) -> list[str]:
    """Extract source-provided antibody labels such as Fab35 from titles."""
    return sorted(
        set(re.findall(r"\b(?:mAb|Fab|Nb)[A-Za-z0-9_-]+\b", text(value), re.IGNORECASE)),
        key=str.casefold,
    )


def _positioned_values(value: object) -> list[str]:
    """Split a positional antigen column while preserving empty positions."""
    raw = text(value)
    if not raw:
        return []
    return [
        "" if part.strip().casefold() in NA else part.strip()
        for part in re.split(r"[;|]", raw)
    ]


def _sabdab_targets(row: dict) -> list[dict]:
    """Return SAbDab2 antigen names with any explicit stable IDs."""
    names = split_values(row.get("antigen_name"))
    types = split_values(row.get("antigen_type"))
    identifier_fields = (
        ("antigen_uniprot", "UniProt"),
        ("antigen_uniprot_id", "UniProt"),
        ("antigen_uniprot_accession", "UniProt"),
        ("antigen_hgnc", "HGNC"),
        ("antigen_hgnc_id", "HGNC"),
        ("antigen_ncbi_gene", "NCBI Gene"),
        ("antigen_ncbi_gene_id", "NCBI Gene"),
        ("antigen_id", ""),
        ("antigen_ids", ""),
        ("uniprot", "UniProt"),
        ("target_external_id", ""),
        ("antigen_identifier", ""),
    )
    positioned_identifiers = []
    for field, namespace in identifier_fields:
        values = _positioned_values(row.get(field))
        if values:
            positioned_identifiers.append(
                (
                    namespace,
                    [
                        value if ":" in value or not namespace else f"{namespace}:{value}"
                        for value in values
                    ],
                )
            )
    if not names:
        for _, values in positioned_identifiers:
            names = [value for value in values if value]
            if names:
                break
    organism = text(
        row.get("antigen_organism")
        or row.get("antigen_species")
        or row.get("agtaxid")
    )
    output = []
    for index in range(len(names)):
        name = names[index]
        antigen_type = types[index].upper() if index < len(types) else ""
        external_id = ""
        for _, values in positioned_identifiers:
            # A singleton identifier is safe only for a singleton antigen.
            # Multi-antigen rows require positional coverage from the source.
            if len(names) == 1 and len(values) == 1:
                external_id = values[0]
                break
            if len(values) == len(names) and index < len(values) and values[index]:
                external_id = values[index]
                break
        # SAbDab2 marks absent antigen fields as NA/NONE.  Other antigen
        # types are retained: the source explicitly records a contact, while
        # target normalization decides whether a name is a curated entity.
        if name and antigen_type not in {"NA", "NONE"}:
            output.append(
                {
                    "name": name,
                    "external_id": external_id,
                    "organism": organism,
                }
            )
    return output


def sabdab2(path: Path, source_url: str = ""):
    """Ingest a SAbDab2 API snapshot plus its bulk instance summary.

    The build downloader stores canonical API records and the all-summary CSV
    together as JSON.  Canonical records provide variable-domain sequences;
    summary rows provide explicit antibody-instance/antigen pairings and PDB
    provenance.  No target is inferred from a structure title or keyword.
    """
    payload = _open_json(path)
    records = payload.get("records", []) if isinstance(payload, dict) else payload
    summary = payload.get("summary", []) if isinstance(payload, dict) else []
    if not isinstance(records, list) or not isinstance(summary, list):
        raise ValueError("SAbDab2 snapshot must contain records and summary arrays")

    summary_by_id: dict[str, list[dict]] = defaultdict(list)
    for row in summary:
        sabdab_id = text(row.get("SABDAB_ID"))
        if sabdab_id:
            summary_by_id[sabdab_id].append(row)

    for record in records:
        sabdab_id = text(record.get("id"))
        if not sabdab_id:
            continue
        heavy_segment = record.get("heavy_like_segment") or {}
        light_segment = record.get("light_like_segment") or {}
        heavy = _sabdab_segment_sequence(heavy_segment)
        light = _sabdab_segment_sequence(light_segment)
        rows = summary_by_id.get(sabdab_id, [])
        if not rows:
            rows = [
                {
                    "INSTANCE": instance.get("id", sabdab_id),
                    "PDB": instance.get("id", ""),
                    "SABDAB_ID": sabdab_id,
                    "organism": "",
                    "type": instance.get("type", ""),
                }
                for instance in record.get("instances", [])
            ] or [{"INSTANCE": sabdab_id, "SABDAB_ID": sabdab_id}]

        structures = sorted(
            {
                _sabdab_pdb_id(row.get("PDB"))
                for row in rows
                if _sabdab_pdb_id(row.get("PDB"))
            }
        )
        first_compound = text(rows[0].get("compound"))
        default_labels = _sabdab_labels(first_compound)
        default_name = default_labels[0] if default_labels else first_compound or sabdab_id
        for row in rows:
            instance_id = text(row.get("INSTANCE")) or sabdab_id
            pdb_id = _sabdab_pdb_id(row.get("PDB"))
            compound = text(row.get("compound")) or default_name
            labels = _sabdab_labels(compound)
            ab = _obs(
                "sabdab2",
                instance_id,
                labels[0] if labels else default_name,
                heavy,
                light,
                cdrh1=heavy_segment.get("cdr1", ""),
                cdrh2=heavy_segment.get("cdr2", ""),
                cdrh3=heavy_segment.get("cdr3", ""),
                cdrl1=light_segment.get("cdr1", ""),
                cdrl2=light_segment.get("cdr2", ""),
                cdrl3=light_segment.get("cdr3", ""),
                organism=text(row.get("organism")),
                format=text(row.get("type")),
                heavy_v=text(heavy_segment.get("v_gene")),
                heavy_j=text(heavy_segment.get("j_gene")),
                light_v=text(light_segment.get("v_gene")),
                light_j=text(light_segment.get("j_gene")),
                structures=structures,
                aliases=[sabdab_id, *labels],
                reference=text(row.get("authors")) or text(row.get("compound")),
                source_url=source_url,
                record_url=(
                    f"https://sabdab.opig.stats.ox.ac.uk/pdb/{pdb_id.lower()}/"
                    if pdb_id
                    else ""
                ),
                link_scope="record" if pdb_id else "source_homepage",
                metadata={
                    "sabdab_id": sabdab_id,
                    "instance_id": instance_id,
                    "pdb_id": pdb_id,
                    "antigen_chain": text(row.get("antigen_chain")),
                    "antigen_type": text(row.get("antigen_type")),
                },
            )
            interactions = [
                InteractionObservation(
                    antibody_id=ab.identity(),
                    source="sabdab2",
                    source_record_id=instance_id,
                    target_raw=target["name"],
                    relationship="binds",
                    evidence="STRUCTURE",
                    reference=text(row.get("compound")) or text(row.get("authors")),
                    record_url=ab.record_url,
                    link_scope=ab.link_scope,
                    note="Explicit antibody-antigen instance pairing in SAbDab2.",
                    target_external_id=target["external_id"],
                    target_organism=target["organism"],
                )
                for target in _sabdab_targets(row)
            ]
            yield ab, interactions


def therasabdab(path: Path, source_url: str = ""):
    reader = _require(
        _open_csv(path), ["Therapeutic", "HeavySequence", "LightSequence", "Target"], "Thera-SAbDab"
    )
    for row in reader:
        name = text(row.get("Therapeutic"))
        if not name:
            continue
        multispecific = _thera_multispecific(row)
        structure_tiers = _thera_structure_tiers(row)
        target_groups = split_values(row.get("Target"))
        arm_specs, arm_mapping_status = _thera_arm_mapping(row, target_groups)
        status = "; ".join(
            x
            for x in [text(row.get("Est. Status")), text(row.get("Highest_Clin_Trial (Feb '25)"))]
            if x
        )
        aliases = split_values(row.get("Alternative Therapeutic Names"))
        metadata = {
            "companies": text(row.get("Companies")),
            "year_proposed": text(row.get("Year Proposed")),
            "genetics": text(row.get("Genetics (Bispecifics delimited with semicolon)")),
            "conditions_approved": text(row.get("Conditions Approved")),
        }
        if multispecific:
            metadata.update(
                {
                    "multispecific": True,
                    "target_assignment_status": arm_mapping_status,
                    "construct_targets": target_groups,
                    "construct_target_context": (
                        "Targets are reported for the construct; arm mapping is unavailable."
                        if arm_mapping_status == "unavailable_no_arm_mapping"
                        else "Targets are explicitly mapped to source-designated sequence arms."
                    ),
                    "quarantine_reason": (
                        "Construct-level targets cannot be assigned to sequence arms without "
                        "an explicit source arm-target mapping."
                    ),
                }
            )
        construct_id = "construct_" + digest("therasabdab", name)
        construct = {
            "id": construct_id,
            "name": name,
            "format": text(row.get("Format")),
            "target_raw": target_groups,
            "target_assignment_status": arm_mapping_status if multispecific else "construct_level",
            "therapeutic_status": status,
            "structure_tiers": structure_tiers,
            "source": "therasabdab",
            "source_record_id": name,
        }

        if multispecific:
            # Emit one sequence observation per available arm.  This keeps
            # arm sequences searchable without presenting the construct as a
            # single antibody with every construct target.
            if not arm_specs:
                arm_specs = [{"role": "unavailable", "heavy": "", "light": "", "label": "arm_1"}]
            construct["arms"] = []
            for arm_spec in arm_specs:
                construct["arms"].append(
                    {
                        "id": f"{construct_id}:{arm_spec['label']}",
                        "label": arm_spec["label"],
                        "role": arm_spec["role"],
                        "has_heavy": bool(arm_spec["heavy"]),
                        "has_light": bool(arm_spec["light"]),
                        "target_assignment_status": arm_mapping_status,
                    }
                )
            for index, arm_spec in enumerate(arm_specs, 1):
                arm_id = f"{construct_id}:{arm_spec['label']}"
                arm = {
                    "id": arm_id,
                    "label": arm_spec["label"],
                    "role": arm_spec["role"],
                    "target_assignment_status": arm_mapping_status,
                }
                if arm_mapping_status == "arm_mapped_explicit_source_field":
                    arm["target_raw"] = arm_spec["target_raw"]
                arm_metadata = dict(metadata)
                arm_metadata["arm_id"] = arm_id
                arm_metadata["arm_label"] = arm_spec["label"]
                arm_name = name if index == 1 else f"{name} [{arm_spec['label']}]"
                arm_aliases = aliases if index == 1 else [name, *aliases]
                ab = _obs(
                    "therasabdab",
                    name,
                    arm_name,
                    arm_spec["heavy"],
                    arm_spec["light"],
                    sequence_accessions=_source_accessions(row),
                    format=text(row.get("Format")),
                    # Therapeutic status and structure matches belong to the
                    # construct unless the source explicitly annotates arms.
                    therapeutic_status="",
                    aliases=arm_aliases,
                    structures=[],
                    structure_tiers={},
                    construct=construct,
                    arm=arm,
                    source_url=source_url,
                    record_url=_thera_record_url(name),
                    link_scope="record",
                    metadata=arm_metadata,
                )
                interactions = []
                if arm_mapping_status == "arm_mapped_explicit_source_field":
                    interactions.append(
                        InteractionObservation(
                            antibody_id=ab.identity(),
                            source="therasabdab",
                            source_record_id=name,
                            target_raw=arm_spec["target_raw"],
                            relationship="targets",
                            evidence="CURATED",
                            reference="Thera-SAbDab therapeutic target annotation",
                            note="Arm target mapped by an explicit source arm-target field.",
                        )
                    )
                yield ab, interactions
            continue

        ab = _obs(
            "therasabdab",
            name,
            name,
            row.get("HeavySequence"),
            row.get("LightSequence"),
            sequence_accessions=_source_accessions(row),
            format=text(row.get("Format")),
            therapeutic_status=status,
            aliases=aliases,
            structures=structure_tiers["100%"],
            structure_tiers=structure_tiers,
            source_url=source_url,
            record_url=_thera_record_url(name),
            link_scope="record",
            metadata=metadata,
        )
        aid = ab.identity()
        interactions = [
            InteractionObservation(
                antibody_id=aid,
                source="therasabdab",
                source_record_id=ab.record_id,
                target_raw=target_group,
                relationship="targets",
                evidence="CURATED",
                reference="Thera-SAbDab therapeutic target annotation",
            )
            for target_group in target_groups
        ]
        yield ab, interactions


def _cov_target(virus: str, epitope: str) -> list[tuple[str, str]]:
    virus = text(virus).replace("_", " ")
    epitope = text(epitope)
    out = [(virus, "source")]
    u = virus.upper()
    if "SARS-COV2" in u or "SARS-COV-2" in u:
        # CoV-AbDab uses free-text epitope annotations.  Parse RBD as a
        # token, and remove explicitly negated occurrences before deciding
        # whether a positive RBD assertion exists.  A substring check would
        # incorrectly turn ``non-RBD`` into a positive RBD target.
        annotation = re.sub(r"[^A-Z0-9]+", " ", epitope.upper()).strip()
        non_rbd = bool(re.search(r"\bNON\s+RBD\b", annotation))
        without_negated_rbd = re.sub(r"\bNON\s+RBD\b", " ", annotation)
        positive_rbd = bool(re.search(r"\bRBD\b", without_negated_rbd))
        mixed_rbd = non_rbd and positive_rbd
        ntd = bool(re.search(r"\bNTD\b", annotation))
        broad_spike = bool(
            re.search(r"\bS\b|\bS1\b|\bSPIKE\b|\bNTD\b", annotation)
            or positive_rbd
            or mixed_rbd
        )
        if positive_rbd and not non_rbd:
            out.append(("SARS-CoV-2 RBD", "source_epitope"))
        elif ntd:
            out.append(("SARS-CoV-2 NTD", "source_epitope"))
        if broad_spike:
            out.append(
                (
                    "SARS-CoV-2 Spike",
                    "ambiguous_epitope" if mixed_rbd else "derived_hierarchy",
                )
            )
    return list(dict.fromkeys(out))


def covabdab(path: Path, source_url: str = ""):
    reader = _require(_open_csv(path), ["Name", "Binds to", "VHorVHH", "VL"], "CoV-AbDab")
    for row in reader:
        name = text(row.get("Name"))
        if not name:
            continue
        ab = _obs(
            "covabdab",
            name,
            name,
            row.get("VHorVHH"),
            row.get("VL"),
            sequence_accessions=_source_accessions(row),
            cdrh3=row.get("CDRH3"),
            cdrl3=row.get("CDRL3"),
            organism=text(row.get("Origin")),
            format=text(row.get("Ab or Nb")),
            heavy_v=text(row.get("Heavy V Gene")),
            heavy_j=text(row.get("Heavy J Gene")),
            light_v=text(row.get("Light V Gene")),
            light_j=text(row.get("Light J Gene")),
            structures=_structures(row.get("Structures", "")),
            reference=text(row.get("Sources")),
            source_url=source_url,
            metadata={
                "date_added": text(row.get("Date Added")),
                "last_updated": text(row.get("Last Updated")),
            },
        )
        aid = ab.identity()
        epi = text(row.get("Protein + Epitope"))
        interactions = []
        relation_fields = [
            ("Binds to", "binds"),
            ("Doesn't Bind to", "does_not_bind"),
            ("Neutralising Vs", "neutralizes"),
            ("Not Neutralising Vs", "does_not_neutralize"),
        ]
        for field, rel in relation_fields:
            for raw in split_values(row.get(field)):
                # Strip weak/parenthetical assay qualifier from entity, retain raw in note.
                entity = re.sub(r"\s*\([^)]*\)\s*$", "", raw).strip()
                targets = (
                    _cov_target(entity, epi)
                    if "bind" in rel
                    else [(entity.replace("_", " "), "source")]
                )
                for t, assertion_origin in targets:
                    interactions.append(
                        InteractionObservation(
                            antibody_id=aid,
                            source="covabdab",
                            source_record_id=ab.record_id,
                            target_raw=t,
                            relationship=rel,
                            evidence="CURATED",
                            reference=ab.reference,
                            epitope=epi,
                            note=raw if raw != entity else "",
                            assertion_origin=assertion_origin,
                        )
                    )
        yield ab, interactions


def _pox_entity(raw: str) -> tuple[str, str]:
    clean = re.sub(r"\([^)]*\)", "", text(raw)).strip()
    parts = [p for p in clean.split("_") if p]
    if not parts:
        return clean, ""
    virus_map = {
        "VACV": "Vaccinia virus",
        "MPXV": "Mpox virus",
        "VARV": "Variola virus",
        "CPXV": "Cowpox virus",
    }
    virus = virus_map.get(parts[0].upper(), parts[0])
    # Retain every source dimension after the virus token. Broadening a
    # strain/location-specific observation to the virus or protein alone is
    # unsafe for retrieval.
    detail = " ".join(parts[1:])
    return f"{virus} {detail}".strip(), clean


def poxabdab(path: Path, source_url: str = ""):
    reader = _require(
        _open_csv(path),
        ["Name", "Binds to(virus[-clade]_strain_location_protein)", "VH or VHH", "VL"],
        "Pox-AbDab",
    )
    for row in reader:
        name = text(row.get("Name"))
        if not name:
            continue
        ab = _obs(
            "poxabdab",
            name,
            name,
            row.get("VH or VHH"),
            row.get("VL"),
            sequence_accessions=_source_accessions(row),
            cdrh3=row.get("CDRH3"),
            cdrl3=row.get("CDRL3"),
            organism=text(row.get("Species origin")),
            format=text(row.get("Construct")),
            heavy_v=text(row.get("Heavy V gene")),
            heavy_j=text(row.get("Heavy J gene")),
            light_v=text(row.get("Light V gene")),
            light_j=text(row.get("Light J gene")),
            structures=_structures(row.get("Solved structures", "")),
            reference=text(row.get("Sources")),
            source_url=source_url,
            metadata={
                "added": text(row.get("Added")),
                "modified": text(row.get("Modified")),
                "epitope_type": text(row.get("Epitope type")),
            },
        )
        aid = ab.identity()
        interactions = []
        relation_fields = [
            ("Binds to(virus[-clade]_strain_location_protein)", "binds"),
            ("Does not bind to(virus[-clade]_strain_location_protein)", "does_not_bind"),
            ("Neutralises(virus[-clade]_strain_location)", "neutralizes"),
            ("Does not neutralise(virus[-clade]_strain_location)", "does_not_neutralize"),
            ("Protects against(virus[-clade]_strain)", "protects"),
            ("Does not protect against(virus[-clade]_strain)", "does_not_protect"),
        ]
        for field, rel in relation_fields:
            for raw in split_values(row.get(field)):
                target, note = _pox_entity(raw)
                interactions.append(
                    InteractionObservation(
                        antibody_id=aid,
                        source="poxabdab",
                        source_record_id=ab.record_id,
                        target_raw=target,
                        relationship=rel,
                        evidence="CURATED",
                        reference=ab.reference,
                        epitope=text(row.get("Epitope type")),
                        note=note,
                    )
                )
        yield ab, interactions


def _iedb_chain_role(chain_type: object) -> str:
    value = text(chain_type).casefold()
    if value in {"h", "heavy"} or "heavy" in value:
        return "VH"
    if value in {"l", "light", "kappa", "lambda"} or any(
        token in value for token in ("light", "kappa", "lambda")
    ):
        return "VL"
    return ""


def _iedb_chain(row: dict, prefix: str) -> tuple[str, dict]:
    role = _iedb_chain_role(row.get(f"{prefix}__type"))
    if not role:
        return "", {}
    regions = {}
    for number in (1, 2, 3):
        curated = text(row.get(f"{prefix}__cdr{number}_curated"))
        calculated = text(row.get(f"{prefix}__cdr{number}_calculated"))
        regions[f"CDR{number}"] = {
            "chosen": calculated or curated,
            "calculated": calculated,
            "curated": curated,
            "calculated_start": text(row.get(f"{prefix}__cdr{number}_start_calculated")),
            "calculated_end": text(row.get(f"{prefix}__cdr{number}_end_calculated")),
            "curated_start": text(row.get(f"{prefix}__cdr{number}_start_curated")),
            "curated_end": text(row.get(f"{prefix}__cdr{number}_end_curated")),
        }
    genes = {}
    for gene in ("v", "d", "j"):
        genes[gene.upper()] = {
            "chosen": text(row.get(f"{prefix}__calculated_{gene}_gene"))
            or text(row.get(f"{prefix}__curated_{gene}_gene")),
            "calculated": text(row.get(f"{prefix}__calculated_{gene}_gene")),
            "curated": text(row.get(f"{prefix}__curated_{gene}_gene")),
        }
    return role, {
        "chain_type": text(row.get(f"{prefix}__type")),
        "organism_iri": text(row.get(f"{prefix}__organism_iri")),
        "protein_iri": text(row.get(f"{prefix}__protein_iri")),
        "source_full_protein": text(row.get(f"{prefix}__protein_sequence")),
        "calculated_v_domain": text(row.get(f"{prefix}__v_domain_calculated")),
        "source_nucleotide": text(row.get(f"{prefix}__nucleotide_sequence")),
        "numbering_scheme": "IMGT" if any(r["calculated"] for r in regions.values()) else "",
        "annotation_origin": "IEDB calculated"
        if any(r["calculated"] for r in regions.values())
        else "IEDB curated",
        "regions": regions,
        "genes": genes,
    }


def _postgres_array(value: object) -> list[str]:
    raw = text(value)
    if not raw or raw == "{}":
        return []
    return [item.strip().strip('"') for item in raw.strip("{}").split(",") if item.strip()]


def _iedb_metric(response: object) -> str:
    value = text(response).upper().replace(" ", "")
    for metric in ("IC50", "EC50", "KOFF", "KON", "KD"):
        if metric in value:
            return metric
    return ""


def _iedb_direct_binding_assay(method: object, response: object, metric: str) -> bool:
    """Return whether an IEDB assay explicitly supports direct binding.

    A positive B-cell result is not, by itself, a molecular binding claim.
    Kinetic/affinity metrics and recognizable binding assay methods are the
    narrow allowlist used for direct ``binds``/``does_not_bind`` edges.
    """
    if metric in {"KD", "KON", "KOFF"}:
        return True
    descriptor = f"{text(method)} {text(response)}".casefold()
    return bool(
        re.search(
            r"(?:surface\s+plasmon|\bspr\b|biolayer\s+interferometry|\bbli\b|"
            r"binding|\belisa\b)",
            descriptor,
        )
    )


def _iedb_support(path: Path) -> tuple[dict[str, list[str]], dict[str, list[dict]]]:
    search_path = path.with_name("iedb-bcr-search.csv")
    measurement_path = path.with_name("iedb-bcell-export.csv")
    if not search_path.exists() or not measurement_path.exists():
        return {}, {}
    assays_by_group: dict[str, set[str]] = {}
    search_rows = _require(
        _open_csv(search_path), ["receptor_group_id", "iedb_assay_ids"], "IEDB bcr_search"
    )
    for row in search_rows:
        group_id = text(row.get("receptor_group_id"))
        if group_id:
            assays_by_group.setdefault(group_id, set()).update(
                _postgres_array(row.get("iedb_assay_ids"))
            )
    rows_by_assay: dict[str, list[dict]] = {}
    measurement_rows = _require(
        _open_csv(measurement_path),
        ["assay_id", "assay__method"],
        "IEDB bcell_export",
    )
    for row in measurement_rows:
        assay_id = text(row.get("assay_id"))
        if assay_id:
            rows_by_assay.setdefault(assay_id, []).append(row)
    return {group: sorted(ids) for group, ids in assays_by_group.items()}, rows_by_assay


def _iedb_interaction(ab, row: dict, group_id: str) -> InteractionObservation | None:
    target_fields = (
        ("epitope__molecule_parent", "epitope__molecule_parent_iri"),
        ("epitope__source_molecule", "epitope__source_molecule_iri"),
    )
    target, target_external_id = next(
        (
            (text(row.get(name)), text(row.get(iri)))
            for name, iri in target_fields
            if text(row.get(name))
        ),
        ("", ""),
    )
    if not target or not target_external_id:
        return None
    response = text(row.get("assay__response_measured"))
    method = text(row.get("assay__method"))
    qualitative = text(row.get("assay__qualitative_measure")).casefold()
    functional = "neutral" in f"{response} {method}".casefold()
    negative = qualitative.startswith("negative")
    positive = qualitative.startswith("positive")
    metric = _iedb_metric(response)
    raw_measurement = text(row.get("assay__quantitative_measurement"))
    if not negative and not positive and not (metric and raw_measurement):
        return None
    if functional:
        relationship = (
            "does_not_neutralize"
            if negative
            else "neutralizes"
        )
    elif not _iedb_direct_binding_assay(method, response, metric):
        return None
    else:
        relationship = "does_not_bind" if negative else "binds"
    assay_id = text(row.get("assay_id"))
    measurements = []
    if metric and raw_measurement:
        measurements.append(
            measurement(
                metric,
                raw_measurement,
                row.get("assay__units"),
                qualifier=text(row.get("assay__measurement_inequality")),
                assay_id=assay_id,
                method=method,
                response_measured=response,
                source_field="assay__quantitative_measurement",
            )
        )
    note = text(row.get("assay__comments"))
    if raw_measurement and not metric:
        note = "; ".join(
            part
            for part in [
                note,
                f"Unmapped quantitative response: {response} {raw_measurement} {text(row.get('assay__units'))}",
            ]
            if part
        )
    return InteractionObservation(
        antibody_id=ab.identity(),
        source="iedb",
        source_record_id=group_id,
        target_raw=target,
        relationship=relationship,
        evidence="MEASURED" if measurements else "CURATED",
        reference=text(row.get("reference__iedb_iri")),
        epitope=text(row.get("epitope__name")),
        assay=" · ".join(part for part in [method, response] if part),
        note=note,
        record_url=text(row.get("assay_id__iedb_iri")),
        link_scope="record",
        measurements=measurements,
        target_external_id=target_external_id,
        assay_ids=[assay_id] if assay_id else [],
        receptor_group_id=group_id,
    )


def iedb(path: Path, source_url: str = ""):
    assays_by_group, rows_by_assay = _iedb_support(path)
    emitted: set[tuple[str, str, str, str, str]] = set()
    reader = _require(
        _open_csv(path),
        [
            "receptor_group_id",
            "receptor__type",
            "epitope__name",
            "chain_1__type",
            "chain_2__type",
        ],
        "IEDB BCR export",
    )
    for row in reader:
        group_id = text(row.get("receptor_group_id"))
        if not group_id:
            continue
        chains = {}
        for prefix in ("chain_1", "chain_2"):
            role, annotation = _iedb_chain(row, prefix)
            if role and role not in chains:
                chains[role] = annotation

        heavy = chains.get("VH", {})
        light = chains.get("VL", {})
        record_url = text(row.get("receptor__group_iri"))
        reference = text(row.get("reference__iedb_iri"))
        name = text(row.get("receptor__reference_name")) or f"IEDB BCR {group_id}"
        nt_provenance = {}
        for role, annotation in chains.items():
            if annotation.get("source_nucleotide"):
                nt_provenance[role] = {
                    "source": "iedb",
                    "source_record_id": group_id,
                    "source_field": next(
                        f"chain_{index}__nucleotide_sequence"
                        for index in (1, 2)
                        if _iedb_chain_role(row.get(f"chain_{index}__type")) == role
                    ),
                    "scope": "full_length_bcr",
                }

        ab = _obs(
            "iedb",
            group_id,
            name,
            heavy.get("calculated_v_domain", ""),
            light.get("calculated_v_domain", ""),
            cdrh1=heavy.get("regions", {}).get("CDR1", {}).get("chosen", ""),
            cdrh2=heavy.get("regions", {}).get("CDR2", {}).get("chosen", ""),
            cdrh3=heavy.get("regions", {}).get("CDR3", {}).get("chosen", ""),
            cdrl1=light.get("regions", {}).get("CDR1", {}).get("chosen", ""),
            cdrl2=light.get("regions", {}).get("CDR2", {}).get("chosen", ""),
            cdrl3=light.get("regions", {}).get("CDR3", {}).get("chosen", ""),
            vh_nt_source=heavy.get("source_nucleotide", ""),
            vl_nt_source=light.get("source_nucleotide", ""),
            sequence_accessions=_source_accessions(row),
            nucleotide_provenance=nt_provenance,
            format=f"IEDB {text(row.get('receptor__type'))}".strip(),
            heavy_v=heavy.get("genes", {}).get("V", {}).get("chosen", ""),
            heavy_d=heavy.get("genes", {}).get("D", {}).get("chosen", ""),
            heavy_j=heavy.get("genes", {}).get("J", {}).get("chosen", ""),
            light_v=light.get("genes", {}).get("V", {}).get("chosen", ""),
            light_j=light.get("genes", {}).get("J", {}).get("chosen", ""),
            chain_annotations=chains,
            reference=reference,
            source_url=source_url,
            record_url=record_url,
            link_scope="record" if record_url else "source_homepage",
            metadata={
                "iedb_receptor_id": text(row.get("receptor__iedb_receptor_id")),
                "epitope_iri": text(row.get("epitope__iedb_iri")),
                "epitope_name": text(row.get("epitope__name")),
                "epitope_source_organism": text(row.get("epitope__source_organism")),
            },
        )
        interactions = []
        support_rows = [
            support
            for assay_id in assays_by_group.get(group_id, [])
            for support in rows_by_assay.get(assay_id, [])
        ]
        for support in support_rows:
            interaction = _iedb_interaction(ab, support, group_id)
            if interaction:
                token = (
                    ab.identity(),
                    group_id,
                    text(support.get("assay_id")),
                    interaction.target_external_id,
                    json.dumps(interaction.measurements, sort_keys=True),
                )
                if token in emitted:
                    continue
                interactions.append(interaction)
                emitted.add(token)
        # A textual source-molecule value is not assay evidence.  In
        # particular, do not manufacture a positive binds edge when the
        # explicit IEDB support exports have no rows for this receptor group.
        yield ab, interactions


def _asd_rows(path: Path) -> list[dict]:
    """Read an authorized ASD export in CSV or JSON form.

    The published ASD schema is tabular, but local users may save a pandas
    export as JSON.  Supporting both forms keeps this adapter useful without
    adding a downloader for the restricted upstream dataset.
    """
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    decoded = raw.decode("utf-8-sig", errors="replace")
    stripped = decoded.lstrip()
    if stripped.startswith("{") or stripped.startswith("["):
        payload = json.loads(decoded)
        if isinstance(payload, dict):
            payload = payload.get("records", payload.get("data", []))
        if not isinstance(payload, list) or not all(isinstance(row, dict) for row in payload):
            raise ValueError("ASD JSON export must contain a list of record objects")
        return payload
    return list(csv.DictReader(io.StringIO(decoded, newline="")))


def _oas_rows(path: Path) -> tuple[dict, list[dict]]:
    """Read one OAS data-unit, whose first line is JSON metadata."""
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    lines = raw.decode("utf-8-sig", errors="replace").splitlines()
    if not lines:
        raise ValueError("OAS data-unit is empty")
    try:
        metadata = json.loads(lines[0])
    except (TypeError, ValueError) as exc:
        raise ValueError("OAS data-unit first line must be JSON metadata") from exc
    if not isinstance(metadata, dict):
        raise ValueError("OAS data-unit metadata must be a JSON object")
    rows = list(csv.DictReader(io.StringIO("\n".join(lines[1:]), newline="")))
    return metadata, rows


def _oas_aa(value: object) -> str:
    """Remove alignment-only gaps before applying the normal AA contract."""
    return text(value).replace("-", "").replace(".", "")


def _oas_sequence(row: dict, *fields: str) -> str:
    for field in fields:
        value = _oas_aa(row.get(field))
        if value:
            try:
                return sequence(value)
            except SequenceNormalizationError:
                return ""
    return ""


def oas(path: Path, source_url: str = ""):
    """Import one explicitly selected OAS data-unit without target inference.

    OAS is a repertoire sequence resource, not an antigen-specificity source.
    Its records therefore contribute searchable sequence observations only;
    publication links, donor/study metadata, and source nucleotide strings are
    preserved, but no interaction is emitted.
    """
    unit_metadata, rows = _oas_rows(path)
    locus = text(unit_metadata.get("Chain") or unit_metadata.get("chain"))
    study = text(unit_metadata.get("Study") or unit_metadata.get("study"))
    author = text(unit_metadata.get("Author") or unit_metadata.get("author"))
    publication_url = text(unit_metadata.get("Link") or unit_metadata.get("link")) or source_url
    species = text(unit_metadata.get("Species") or unit_metadata.get("species"))
    for index, row in enumerate(rows, 1):
        row_locus = text(row.get("locus")) or locus
        heavy = _oas_sequence(
            row,
            "heavy_sequence_alignment_aa",
            "heavy_v_sequence_alignment_aa",
        )
        light = _oas_sequence(
            row,
            "light_sequence_alignment_aa",
            "light_v_sequence_alignment_aa",
        )
        if not heavy and not light:
            if row_locus.upper().startswith("IGH"):
                heavy = _oas_sequence(row, "sequence_alignment_aa", "v_sequence_alignment_aa")
            elif row_locus.upper().startswith(("IGK", "IGL")):
                light = _oas_sequence(row, "sequence_alignment_aa", "v_sequence_alignment_aa")
        if not heavy and not light:
            continue
        vh_nt = ""
        vl_nt = ""
        raw_nt = text(row.get("sequence"))
        if raw_nt:
            try:
                if row_locus.upper().startswith("IGH"):
                    vh_nt = nucleotide_sequence(raw_nt)
                elif row_locus.upper().startswith(("IGK", "IGL")):
                    vl_nt = nucleotide_sequence(raw_nt)
            except SequenceNormalizationError:
                pass
        sequence_id = text(row.get("sequence_id")) or text(row.get("sequence"))
        record_id = f"{study or 'oas'}:{sequence_id or index}"
        row_metadata = {
            "oas_study": study,
            "oas_author": author,
            "oas_unit_metadata": unit_metadata,
            "oas_locus": row_locus,
            "oas_isotype": text(row.get("c_region")),
            "oas_species": species,
            "oas_sequence_only": True,
        }
        row_metadata = {key: value for key, value in row_metadata.items() if value not in ("", {})}
        nt_provenance = {}
        if vh_nt:
            nt_provenance["VH"] = {
                "source": "oas",
                "source_record_id": record_id,
                "scope": "full_length_repertoire_sequence",
                "source_field": "sequence",
            }
        if vl_nt:
            nt_provenance["VL"] = {
                "source": "oas",
                "source_record_id": record_id,
                "scope": "full_length_repertoire_sequence",
                "source_field": "sequence",
            }
        ab = _obs(
            "oas",
            record_id,
            f"OAS {record_id}",
            heavy,
            light,
            cdrh1=_oas_sequence(row, "cdr1_aa") if heavy else "",
            cdrh2=_oas_sequence(row, "cdr2_aa") if heavy else "",
            cdrh3=_oas_sequence(row, "cdr3_aa") if heavy else "",
            cdrl1=_oas_sequence(row, "cdr1_aa") if light and not heavy else "",
            cdrl2=_oas_sequence(row, "cdr2_aa") if light and not heavy else "",
            cdrl3=_oas_sequence(row, "cdr3_aa") if light and not heavy else "",
            vh_nt_source=vh_nt,
            vl_nt_source=vl_nt,
            nucleotide_provenance=nt_provenance,
            organism=species,
            format=f"OAS {row_locus}".strip(),
            reference=publication_url,
            source_url=source_url,
            record_url=publication_url,
            link_scope="record" if publication_url else "source_homepage",
            metadata=row_metadata,
        )
        yield ab, []


def _asd_metadata(row: dict) -> dict:
    metadata = {}
    nested = row.get("metadata")
    if isinstance(nested, dict):
        metadata.update(nested)
    elif text(nested):
        try:
            decoded = json.loads(text(nested))
        except (TypeError, ValueError):
            decoded = {}
        if isinstance(decoded, dict):
            metadata.update(decoded)
    for field, value in row.items():
        if field.startswith("metadata."):
            metadata[field.removeprefix("metadata.")] = value
    return metadata


def _asd_value(row: dict, metadata: dict, *names: str) -> object:
    for name in names:
        if name in row and text(row.get(name)):
            return row.get(name)
        if name in metadata and text(metadata.get(name)):
            return metadata.get(name)
    return ""


def _asd_json_value(value: object) -> object:
    if isinstance(value, str) and value.lstrip().startswith(("{", "[")):
        try:
            return json.loads(value)
        except (TypeError, ValueError):
            pass
    return value


def _asd_cdrs(metadata: dict, chain: str) -> dict[str, str]:
    numbering = _asd_json_value(
        _asd_value(
            {},
            metadata,
            f"{chain}_riot_numbering",
            f"{chain}_numbering",
        )
    )
    if isinstance(numbering, dict):
        return {
            f"cdr{chain[0].lower()}{index}": text(
                numbering.get(f"cdr{index}_aa") or numbering.get(f"CDR{index}")
            )
            for index in (1, 2, 3)
        }
    raw = text(numbering)
    output = {}
    for index in (1, 2, 3):
        match = re.search(
            rf"cdr{index}_aa\s*[:=]\s*['\"]?([A-Za-z*.-]+)", raw, re.IGNORECASE
        )
        if match:
            output[f"cdr{chain[0].lower()}{index}"] = match.group(1)
    return output


def _asd_bool(value: object) -> bool:
    return text(value).casefold() in {"1", "true", "yes", "y"}


def _asd_identifier(value: object, namespace: str) -> str:
    raw = text(value)
    if not raw:
        return ""
    if namespace == "UniProt":
        match = re.search(r"(?:UNIPROT(?:KB)?[:/]|uniprot\.org/(?:uniprotkb/)?)([A-Z0-9_-]+)", raw, re.I)
        raw = match.group(1) if match else raw.split()[0]
        return f"UniProt:{raw.upper()}"
    if namespace == "HGNC":
        match = re.search(r"HGNC[:/]?(\d+)", raw, re.I)
        raw = match.group(1) if match else raw
        return f"HGNC:{raw}"
    if namespace == "NCBI Gene":
        match = re.search(r"(?:NCBI[_ ]?GENE|ENTREZ|GENEID)[:/]?(\d+)", raw, re.I)
        raw = match.group(1) if match else raw
        return f"NCBI_GENE:{raw}"
    if namespace == "ChEBI":
        match = re.search(r"CHEBI[:/]?(\d+)", raw, re.I)
        raw = match.group(1) if match else raw
        return f"ChEBI:{raw}"
    match = re.search(r"\b([0-9][A-Za-z0-9]{3})\b", raw)
    return f"PDB:{(match.group(1) if match else raw).upper()}"


def _asd_target(metadata: dict) -> tuple[str, str]:
    for name, namespace in (
        ("target_uniprot", "UniProt"),
        ("target_hgnc", "HGNC"),
        ("target_ncbi_gene", "NCBI Gene"),
        ("target_chebi", "ChEBI"),
        ("target_pdb", "PDB"),
    ):
        identifier = _asd_identifier(_asd_value({}, metadata, name), namespace)
        if identifier:
            return identifier, text(_asd_value({}, metadata, "target_name")) or identifier
    return "", text(_asd_value({}, metadata, "target_name"))


def _asd_affinity(row: dict, metadata: dict) -> tuple[str, list[dict], str]:
    affinity_type = text(_asd_value(row, metadata, "affinity_type", "measurement_type"))
    raw_affinity = text(_asd_value(row, metadata, "affinity", "measurement"))
    if not raw_affinity:
        return "binds", [], affinity_type
    if affinity_type.casefold() in {"bool", "boolean", "binary"}:
        negative = raw_affinity.casefold() in {
            "0",
            "0.0",
            "false",
            "no",
            "negative",
            "non-binding",
        }
        return ("does_not_bind" if negative else "binds"), [], affinity_type

    metric = ""
    for candidate, pattern in (
        ("KD", r"\bk\s*[_ -]?d\b|dissociation"),
        ("KON", r"\bk\s*[_ -]?on\b|association"),
        ("KOFF", r"\bk\s*[_ -]?off\b"),
        ("IC50", r"\bic\s*[_ -]?50\b"),
        ("EC50", r"\bec\s*[_ -]?50\b"),
    ):
        if re.search(pattern, affinity_type, re.IGNORECASE):
            metric = candidate
            break
    if not metric:
        return "binds", [], affinity_type
    match = re.match(
        r"\s*(<=|>=|<|>|~|=)?\s*([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)\s*(.*)",
        raw_affinity,
    )
    if not match:
        return "binds", [], affinity_type
    qualifier, numeric, unit = match.groups()
    return (
        "binds",
        [measurement(metric, numeric, unit, qualifier or "", source_field="affinity")],
        affinity_type,
    )


def naturalantibody_asd(path: Path, source_url: str = ""):
    """Import an authorized NaturalAntibody ASD export for local/private builds.

    ASD is not downloaded by PAIRS and remains disabled in the public config.
    Target IDs are preferred over names, while every source row and upstream
    URL remains independently represented in the normal provenance model.
    """
    rows = _asd_rows(path)
    for index, row in enumerate(rows, 1):
        metadata = _asd_metadata(row)
        dataset = text(_asd_value(row, metadata, "dataset")) or "unknown"
        source_record_id = text(
            _asd_value(row, metadata, "record_id", "antibody_id", "sample_id", "id")
        ) or str(index)
        record_id = f"{dataset}:{source_record_id}"
        heavy = _asd_value(row, metadata, "heavy_sequence", "heavy")
        light = _asd_value(row, metadata, "light_sequence", "light")
        cdrs = {}
        cdrs.update(_asd_cdrs(metadata, "heavy"))
        cdrs.update(_asd_cdrs(metadata, "light"))
        scfv = _asd_bool(_asd_value(row, metadata, "scfv"))
        format_value = text(_asd_value(row, metadata, "format", "antibody_format"))
        if scfv and not format_value:
            format_value = "scFv"
        target_external_id, target_name = _asd_target(metadata)
        antigen_sequence = text(_asd_value(row, metadata, "antigen_sequence"))
        if not target_name and antigen_sequence:
            target_name = "ASD antigen sequence " + digest(antigen_sequence, length=12)
        if not target_name:
            continue
        relationship, measurements, affinity_type = _asd_affinity(row, metadata)
        upstream_url = text(_asd_value(row, metadata, "source_url")) or source_url
        confidence = text(_asd_value(row, metadata, "confidence"))
        observation_metadata = {
            "asd_dataset": dataset,
            "confidence": confidence,
            "affinity_type": affinity_type,
            "target_uniprot": text(_asd_value(row, metadata, "target_uniprot")),
            "target_pdb": text(_asd_value(row, metadata, "target_pdb")),
            "target_name": target_name,
            "antigen_sequence": antigen_sequence,
        }
        observation_metadata = {key: value for key, value in observation_metadata.items() if value}
        name = text(_asd_value(row, metadata, "antibody_name", "name")) or f"ASD {record_id}"
        ab = _obs(
            "naturalantibody_asd",
            record_id,
            name,
            heavy,
            light,
            cdrh1=cdrs.get("cdrh1", ""),
            cdrh2=cdrs.get("cdrh2", ""),
            cdrh3=cdrs.get("cdrh3", ""),
            cdrl1=cdrs.get("cdrl1", ""),
            cdrl2=cdrs.get("cdrl2", ""),
            cdrl3=cdrs.get("cdrl3", ""),
            format=format_value,
            sequence_architecture="scfv" if scfv else "",
            reference=upstream_url,
            source_url=source_url,
            record_url=upstream_url,
            link_scope="record" if upstream_url else "source_homepage",
            aliases=[source_record_id],
            metadata=observation_metadata,
        )
        note = "; ".join(
            part for part in [
                f"ASD confidence: {confidence}" if confidence else "",
                f"ASD dataset: {dataset}",
            ] if part
        )
        yield ab, [
            InteractionObservation(
                antibody_id=ab.identity(),
                source="naturalantibody_asd",
                source_record_id=record_id,
                target_raw=target_name,
                relationship=relationship,
                evidence="MEASURED" if measurements else "CURATED",
                reference=upstream_url,
                assay=affinity_type,
                note=note,
                record_url=upstream_url,
                link_scope="record" if upstream_url else "source_homepage",
                measurements=measurements,
                target_external_id=target_external_id,
                target_organism=text(
                    _asd_value(row, metadata, "target_organism", "source_organism")
                ),
            )
        ]


def _offline_targeted_export(
    path: Path,
    source: str,
    source_url: str,
    id_names: tuple[str, ...],
    antibody_name_names: tuple[str, ...],
    target_name_names: tuple[str, ...],
):
    """Normalize a curated target-linked CSV/JSON export for local builds."""
    for index, row in enumerate(_asd_rows(path), 1):
        metadata = _asd_metadata(row)
        for field in (
            "target_uniprot",
            "target_hgnc",
            "target_ncbi_gene",
            "target_chebi",
            "target_pdb",
            "target_organism",
            "source_url",
        ):
            value = _asd_value(row, metadata, field)
            if value:
                metadata.setdefault(field, value)
        dataset = text(_asd_value(row, metadata, "dataset")) or source
        external_record_id = text(_asd_value(row, metadata, *id_names)) or str(index)
        record_id = f"{dataset}:{external_record_id}"
        metadata["source_record_id"] = external_record_id
        target_name = text(_asd_value(row, metadata, *target_name_names))
        target_external_id, target_from_id = _asd_target(metadata)
        target_name = target_name or target_from_id
        if not target_name:
            continue
        cdrs = {}
        cdrs.update(_asd_cdrs(metadata, "heavy"))
        cdrs.update(_asd_cdrs(metadata, "light"))
        heavy = _asd_value(
            row,
            metadata,
            "heavy_sequence",
            "vh_sequence",
            "heavy_chain_sequence",
            "VH",
            "VH_sequence",
        )
        light = _asd_value(
            row,
            metadata,
            "light_sequence",
            "vl_sequence",
            "light_chain_sequence",
            "VL",
            "VL_sequence",
        )
        format_value = text(
            _asd_value(row, metadata, "format", "antibody_format", "ig_format")
        )
        scfv = _asd_bool(_asd_value(row, metadata, "scfv"))
        if scfv and not format_value:
            format_value = "scFv"
        upstream_url = text(_asd_value(row, metadata, "record_url", "source_url")) or source_url
        name = text(_asd_value(row, metadata, *antibody_name_names)) or f"{source} {record_id}"
        aliases = []
        for field in ("synonyms", "antibody_synonyms", "aliases"):
            aliases.extend(split_values(_asd_value(row, metadata, field), separators=r"[;,|]"))
        structures = _structures(
            _asd_value(row, metadata, "structures", "pdb", "cross_references", "crossrefs")
        )
        observation_metadata = {
            "source_record_id": external_record_id,
            "target_name": target_name,
            "target_uniprot": text(metadata.get("target_uniprot")),
            "target_hgnc": text(metadata.get("target_hgnc")),
            "target_ncbi_gene": text(metadata.get("target_ncbi_gene")),
            "target_chebi": text(metadata.get("target_chebi")),
            "target_pdb": text(metadata.get("target_pdb")),
        }
        observation_metadata = {key: value for key, value in observation_metadata.items() if value}
        ab = _obs(
            source,
            record_id,
            name,
            heavy,
            light,
            cdrh1=cdrs.get("cdrh1", ""),
            cdrh2=cdrs.get("cdrh2", ""),
            cdrh3=cdrs.get("cdrh3", ""),
            cdrl1=cdrs.get("cdrl1", ""),
            cdrl2=cdrs.get("cdrl2", ""),
            cdrl3=cdrs.get("cdrl3", ""),
            format=format_value,
            sequence_architecture="scfv" if scfv else "",
            therapeutic_status=text(
                _asd_value(row, metadata, "therapeutic_status", "status", "development_status")
            ),
            aliases=sorted(set(aliases), key=str.casefold),
            structures=structures,
            reference=text(_asd_value(row, metadata, "reference", "publication", "publications")),
            source_url=source_url,
            record_url=upstream_url,
            link_scope="record" if upstream_url else "source_homepage",
            metadata=observation_metadata,
        )
        yield ab, [
            InteractionObservation(
                antibody_id=ab.identity(),
                source=source,
                source_record_id=record_id,
                target_raw=target_name,
                relationship="binds",
                evidence="STRUCTURE" if structures else "CURATED",
                reference=ab.reference,
                epitope=text(_asd_value(row, metadata, "epitope", "epitope_description")),
                note=f"{source} curated target assignment",
                record_url=upstream_url,
                link_scope="record" if upstream_url else "source_homepage",
                target_external_id=target_external_id,
                target_organism=text(
                    _asd_value(row, metadata, "target_organism", "target_species", "species")
                ),
            )
        ]


def abcd(path: Path, source_url: str = ""):
    """Import an authorized ABCD CSV/JSON export for local/private builds."""
    yield from _offline_targeted_export(
        path,
        "abcd",
        source_url,
        ("abcd_id", "antibody_id", "id"),
        ("antibody_name", "common_name", "name"),
        ("target_name", "target", "antigen_name"),
    )


def imgt_mabdb(path: Path, source_url: str = ""):
    """Import an authorized IMGT/mAb-DB CSV/JSON export locally."""
    yield from _offline_targeted_export(
        path,
        "imgt_mabdb",
        source_url,
        ("imgt_mab_id", "mab_id", "antibody_id", "id"),
        ("inn", "common_name", "proprietary_name", "antibody_name", "name"),
        ("target_name", "specificity", "target", "antigen_name"),
    )


def nanosaurus(path: Path, source_url: str = ""):
    """Import public NanoSaurus nanobodies and their explicit antigen links."""
    payload = _open_json(path)
    rows = payload.get("results", []) if isinstance(payload, dict) else payload
    if not isinstance(rows, list):
        raise ValueError("NanoSaurus snapshot must contain a results array")
    for index, row in enumerate(rows, 1):
        if not isinstance(row, dict):
            continue
        nanobody = row.get("nanobody") or {}
        if not isinstance(nanobody, dict):
            continue
        nanosaurus_id = text(row.get("nanosaurus_id")) or str(index)
        sequence_value = text(nanobody.get("sequence"))
        cdrs = nanobody.get("cdrs") or {}
        if not sequence_value:
            continue
        validation = split_values(nanobody.get("validation_levels"), separators=r"[,;|]")
        structures = _structures(nanobody.get("structures"))
        source_record_url = f"https://www.nanosaurus.org/nanosaurus/{nanosaurus_id}"
        metadata = {
            "nanosaurus_id": nanosaurus_id,
            "protein_library": nanobody.get("protein_library") or {},
            "discovery_strategies": nanobody.get("discovery_strategies") or [],
            "validation_levels": validation,
            "antibody_affinity": text(nanobody.get("antibody_affinity")),
            "pharmacological_property": text(nanobody.get("pharmacological_property")),
            "submitter_orcid_id": text(row.get("submitter_orcid_id")),
        }
        metadata = {key: value for key, value in metadata.items() if value not in ("", [], {})}
        aliases = split_values(nanobody.get("alternative_names"), separators=r"[,;|]")
        ab = _obs(
            "nanosaurus",
            nanosaurus_id,
            text(nanobody.get("common_name")) or nanosaurus_id,
            sequence_value,
            cdrh1=cdrs.get("cdr_1", ""),
            cdrh2=cdrs.get("cdr_2", ""),
            cdrh3=cdrs.get("cdr_3", ""),
            organism=text((nanobody.get("protein_library") or {}).get("antibody_source")),
            format="VHH nanobody",
            sequence_architecture="vhh",
            aliases=aliases,
            structures=structures,
            source_url=source_url or "https://www.nanosaurus.org/",
            record_url=source_record_url,
            link_scope="record",
            metadata=metadata,
        )
        interactions = []
        antigens = row.get("antigens") or []
        if not isinstance(antigens, list):
            antigens = [antigens]
        for antigen in antigens:
            if not isinstance(antigen, dict):
                continue
            target_name = text(antigen.get("antigen_name"))
            target_external_id = ""
            if text(antigen.get("antigen_uniprot_code")):
                target_external_id = _asd_identifier(
                    antigen.get("antigen_uniprot_code"), "UniProt"
                )
            elif text(antigen.get("antigen_pdb_code")):
                target_external_id = _asd_identifier(antigen.get("antigen_pdb_code"), "PDB")
            target_name = target_name or target_external_id
            if not target_name:
                continue
            interactions.append(
                InteractionObservation(
                    antibody_id=ab.identity(),
                    source="nanosaurus",
                    source_record_id=nanosaurus_id,
                    target_raw=target_name,
                    relationship="binds",
                    evidence="STRUCTURE" if structures or any(
                        "structural" in level.casefold() for level in validation
                    ) else "CURATED",
                    note="NanoSaurus public curated nanobody-antigen assignment",
                    record_url=source_record_url,
                    link_scope="record",
                    target_external_id=target_external_id,
                    target_organism=text(antigen.get("antigen_organism")),
                )
            )
        yield ab, interactions


ADAPTERS = {
    "plabdab": plabdab,
    "therasabdab": therasabdab,
    "covabdab": covabdab,
    "poxabdab": poxabdab,
    "iedb": iedb,
    "sabdab2": sabdab2,
    "naturalantibody_asd": naturalantibody_asd,
    "abcd": abcd,
    "imgt_mabdb": imgt_mabdb,
    "oas": oas,
    "nanosaurus": nanosaurus,
}
