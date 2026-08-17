from __future__ import annotations

import csv
import gzip
import io
import re
from urllib.parse import quote
from pathlib import Path

from .model import (
    AntibodyObservation,
    InteractionObservation,
    SequenceNormalizationError,
    digest,
    sequence,
    split_values,
    text,
)

csv.field_size_limit(50_000_000)

THERA_SUMMARY_URL = "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/therasabdab/therasummary/?INN="


def _open_csv(path: Path):
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    decoded = raw.decode("utf-8-sig", errors="replace")
    return csv.DictReader(io.StringIO(decoded, newline=""))


def _require(reader, fields: list[str], source: str):
    available=set(reader.fieldnames or [])
    missing=[f for f in fields if f not in available]
    if missing:
        raise ValueError(f"{source} schema changed; missing columns: {', '.join(missing)}")
    return reader


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

    if any(text(row.get(field)) for field in ("HeavySequence(ifbispec)", "LightSequence(ifbispec)")):
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
    invalid_sequences = {}
    normalized = {}
    for field, value in {
        "heavy": heavy,
        "light": light,
        "cdrh3": kw.pop("cdrh3", ""),
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
    if invalid_sequences:
        metadata.update(
            {
                "sequence_quarantine": True,
                "sequence_quarantine_fields": invalid_sequences,
            }
        )
    return AntibodyObservation(
        source=source,
        record_id=text(rid) or text(name),
        name=text(name) or text(rid) or "Unnamed antibody",
        heavy=normalized["heavy"],
        light=normalized["light"],
        cdrh3=normalized["cdrh3"],
        cdrl3=normalized["cdrl3"],
        metadata=metadata,
        **kw,
    )


def plabdab(path: Path, source_url: str = ""):
    reader=_require(_open_csv(path),["ID","heavy_sequence","light_sequence","targets_mentioned"],"PLAbDab")
    for row in reader:
        target_group = text(row.get("targets_mentioned"))
        if not target_group:
            continue
        ab = _obs(
            "plabdab", row.get("ID"), row.get("ID"), row.get("heavy_sequence"), row.get("light_sequence"),
            organism=text(row.get("organism")),
            format="paired antibody" if text(row.get("light_sequence")) else "antibody",
            reference=text(row.get("reference_title")),
            source_url=source_url,
            metadata={"pairing": text(row.get("pairing")), "update_date": text(row.get("update_date"))},
        )
        aid = ab.identity()
        interactions = []
        for mentioned_target in split_values(target_group, separators=r";"):
            interactions.append(InteractionObservation(
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
            ))
        yield ab, interactions


def therasabdab(path: Path, source_url: str = ""):
    reader=_require(_open_csv(path),["Therapeutic","HeavySequence","LightSequence","Target"],"Thera-SAbDab")
    for row in reader:
        name = text(row.get("Therapeutic"))
        if not name:
            continue
        multispecific = _thera_multispecific(row)
        structure_tiers = _thera_structure_tiers(row)
        target_groups = split_values(row.get("Target"))
        arm_specs, arm_mapping_status = _thera_arm_mapping(row, target_groups)
        status = "; ".join(x for x in [text(row.get("Est. Status")), text(row.get("Highest_Clin_Trial (Feb '25)"))] if x)
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
                    "therasabdab", name, arm_name, arm_spec["heavy"], arm_spec["light"],
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
                    interactions.append(InteractionObservation(
                        antibody_id=ab.identity(),
                        source="therasabdab",
                        source_record_id=name,
                        target_raw=arm_spec["target_raw"],
                        relationship="targets",
                        evidence="CURATED",
                        reference="Thera-SAbDab therapeutic target annotation",
                        note="Arm target mapped by an explicit source arm-target field.",
                    ))
                yield ab, interactions
            continue

        ab = _obs(
            "therasabdab", name, name, row.get("HeavySequence"), row.get("LightSequence"),
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
    out=[(virus, "source")]
    u=virus.upper()
    e=epitope.upper()
    if "SARS-COV2" in u or "SARS-COV-2" in u:
        if "RBD" in e:
            out.append(("SARS-CoV-2 RBD", "source_epitope"))
        elif "NTD" in e:
            out.append(("SARS-CoV-2 NTD", "source_epitope"))
        if re.search(r"(^|[; /])S($|[; /])", e) or "SPIKE" in e or "RBD" in e or "NTD" in e:
            out.append(("SARS-CoV-2 Spike", "derived_hierarchy"))
    return list(dict.fromkeys(out))


def covabdab(path: Path, source_url: str = ""):
    reader=_require(_open_csv(path),["Name","Binds to","VHorVHH","VL"],"CoV-AbDab")
    for row in reader:
        name=text(row.get("Name"))
        if not name: continue
        ab=_obs(
            "covabdab", name, name, row.get("VHorVHH"), row.get("VL"),
            cdrh3=row.get("CDRH3"), cdrl3=row.get("CDRL3"),
            organism=text(row.get("Origin")), format=text(row.get("Ab or Nb")),
            heavy_v=text(row.get("Heavy V Gene")), heavy_j=text(row.get("Heavy J Gene")),
            light_v=text(row.get("Light V Gene")), light_j=text(row.get("Light J Gene")),
            structures=_structures(row.get("Structures", "")), reference=text(row.get("Sources")),
            source_url=source_url,
            metadata={"date_added": text(row.get("Date Added")), "last_updated": text(row.get("Last Updated"))},
        )
        aid=ab.identity(); epi=text(row.get("Protein + Epitope")); interactions=[]
        relation_fields=[
            ("Binds to","binds"),("Doesn't Bind to","does_not_bind"),
            ("Neutralising Vs","neutralizes"),("Not Neutralising Vs","does_not_neutralize")]
        for field,rel in relation_fields:
            for raw in split_values(row.get(field)):
                # Strip weak/parenthetical assay qualifier from entity, retain raw in note.
                entity=re.sub(r"\s*\([^)]*\)\s*$", "", raw).strip()
                targets=(
                    _cov_target(entity, epi)
                    if "bind" in rel
                    else [(entity.replace("_", " "), "source")]
                )
                for t, assertion_origin in targets:
                    interactions.append(InteractionObservation(
                        antibody_id=aid, source="covabdab", source_record_id=ab.record_id,
                        target_raw=t, relationship=rel, evidence="CURATED", reference=ab.reference,
                        epitope=epi, note=raw if raw != entity else "",
                        assertion_origin=assertion_origin,
                    ))
        yield ab, interactions


def _pox_entity(raw: str) -> tuple[str,str]:
    clean=re.sub(r"\([^)]*\)", "", text(raw)).strip()
    parts=[p for p in clean.split("_") if p]
    if not parts: return clean, ""
    virus_map={"VACV":"Vaccinia virus","MPXV":"Mpox virus","VARV":"Variola virus","CPXV":"Cowpox virus"}
    virus=virus_map.get(parts[0].upper(), parts[0])
    # Retain every source dimension after the virus token. Broadening a
    # strain/location-specific observation to the virus or protein alone is
    # unsafe for retrieval.
    detail=" ".join(parts[1:])
    return f"{virus} {detail}".strip(), clean


def poxabdab(path: Path, source_url: str = ""):
    reader=_require(_open_csv(path),["Name","Binds to(virus[-clade]_strain_location_protein)","VH or VHH","VL"],"Pox-AbDab")
    for row in reader:
        name=text(row.get("Name"))
        if not name: continue
        ab=_obs(
            "poxabdab", name, name, row.get("VH or VHH"), row.get("VL"),
            cdrh3=row.get("CDRH3"), cdrl3=row.get("CDRL3"),
            organism=text(row.get("Species origin")), format=text(row.get("Construct")),
            heavy_v=text(row.get("Heavy V gene")), heavy_j=text(row.get("Heavy J gene")),
            light_v=text(row.get("Light V gene")), light_j=text(row.get("Light J gene")),
            structures=_structures(row.get("Solved structures", "")), reference=text(row.get("Sources")),
            source_url=source_url,
            metadata={"added": text(row.get("Added")), "modified": text(row.get("Modified")), "epitope_type": text(row.get("Epitope type"))},
        )
        aid=ab.identity(); interactions=[]
        relation_fields=[
            ("Binds to(virus[-clade]_strain_location_protein)","binds"),
            ("Does not bind to(virus[-clade]_strain_location_protein)","does_not_bind"),
            ("Neutralises(virus[-clade]_strain_location)","neutralizes"),
            ("Does not neutralise(virus[-clade]_strain_location)","does_not_neutralize"),
            ("Protects against(virus[-clade]_strain)","protects"),
            ("Does not protect against(virus[-clade]_strain)","does_not_protect"),
        ]
        for field,rel in relation_fields:
            for raw in split_values(row.get(field)):
                target, note=_pox_entity(raw)
                interactions.append(InteractionObservation(
                    antibody_id=aid, source="poxabdab", source_record_id=ab.record_id,
                    target_raw=target, relationship=rel, evidence="CURATED", reference=ab.reference,
                    epitope=text(row.get("Epitope type")), note=note,
                ))
        yield ab, interactions

ADAPTERS={"plabdab":plabdab,"therasabdab":therasabdab,"covabdab":covabdab,"poxabdab":poxabdab}
