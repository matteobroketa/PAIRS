from __future__ import annotations

import csv
import gzip
import io
import re
from pathlib import Path

from .model import AntibodyObservation, InteractionObservation, sequence, split_values, text

csv.field_size_limit(50_000_000)


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


def _obs(source, rid, name, heavy="", light="", **kw):
    return AntibodyObservation(
        source=source,
        record_id=text(rid) or text(name),
        name=text(name) or text(rid) or "Unnamed antibody",
        heavy=sequence(heavy),
        light=sequence(light),
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
        status = "; ".join(x for x in [text(row.get("Est. Status")), text(row.get("Highest_Clin_Trial (Feb '25)"))] if x)
        aliases = split_values(row.get("Alternative Therapeutic Names"))
        ab = _obs(
            "therasabdab", name, name, row.get("HeavySequence"), row.get("LightSequence"),
            format=text(row.get("Format")),
            therapeutic_status=status,
            aliases=aliases,
            structures=_structures(";".join([
                text(row.get("100% SI Structure")), text(row.get("99% SI Structure")), text(row.get("95-98% SI Structure"))
            ])),
            source_url=source_url,
            metadata={
                "companies": text(row.get("Companies")),
                "year_proposed": text(row.get("Year Proposed")),
                "genetics": text(row.get("Genetics (Bispecifics delimited with semicolon)")),
                "conditions_approved": text(row.get("Conditions Approved")),
            },
        )
        aid = ab.identity()
        interactions=[]
        for target_group in split_values(row.get("Target")):
            interactions.append(InteractionObservation(
                antibody_id=aid,
                source="therasabdab",
                source_record_id=ab.record_id,
                target_raw=target_group,
                relationship="targets",
                evidence="CURATED",
                reference="Thera-SAbDab therapeutic target annotation",
            ))
        yield ab, interactions


def _cov_target(virus: str, epitope: str) -> list[str]:
    virus = text(virus).replace("_", " ")
    epitope = text(epitope)
    out=[virus]
    u=virus.upper()
    e=epitope.upper()
    if "SARS-COV2" in u or "SARS-COV-2" in u:
        if "RBD" in e:
            out.append("SARS-CoV-2 RBD")
        elif "NTD" in e:
            out.append("SARS-CoV-2 NTD")
        if re.search(r"(^|[; /])S($|[; /])", e) or "SPIKE" in e or "RBD" in e or "NTD" in e:
            out.append("SARS-CoV-2 Spike")
    return list(dict.fromkeys(out))


def covabdab(path: Path, source_url: str = ""):
    reader=_require(_open_csv(path),["Name","Binds to","VHorVHH","VL"],"CoV-AbDab")
    for row in reader:
        name=text(row.get("Name"))
        if not name: continue
        ab=_obs(
            "covabdab", name, name, row.get("VHorVHH"), row.get("VL"),
            cdrh3=sequence(row.get("CDRH3")), cdrl3=sequence(row.get("CDRL3")),
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
                targets=_cov_target(entity, epi) if "bind" in rel else [entity.replace("_", " ")]
                for t in targets:
                    interactions.append(InteractionObservation(
                        antibody_id=aid, source="covabdab", source_record_id=ab.record_id,
                        target_raw=t, relationship=rel, evidence="CURATED", reference=ab.reference,
                        epitope=epi, note=raw if raw != entity else "",
                    ))
        yield ab, interactions


def _pox_entity(raw: str) -> tuple[str,str]:
    clean=re.sub(r"\([^)]*\)", "", text(raw)).strip()
    parts=[p for p in clean.split("_") if p]
    if not parts: return clean, ""
    virus_map={"VACV":"Vaccinia virus","MPXV":"Mpox virus","VARV":"Variola virus","CPXV":"Cowpox virus"}
    virus=virus_map.get(parts[0].upper(), parts[0])
    protein=parts[-1] if len(parts)>=2 and re.fullmatch(r"[A-Za-z]+\d+[A-Za-z]*", parts[-1]) else ""
    return (f"{virus} {protein}".strip() if protein else virus), clean


def poxabdab(path: Path, source_url: str = ""):
    reader=_require(_open_csv(path),["Name","Binds to(virus[-clade]_strain_location_protein)","VH or VHH","VL"],"Pox-AbDab")
    for row in reader:
        name=text(row.get("Name"))
        if not name: continue
        ab=_obs(
            "poxabdab", name, name, row.get("VH or VHH"), row.get("VL"),
            cdrh3=sequence(row.get("CDRH3")), cdrl3=sequence(row.get("CDRL3")),
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
