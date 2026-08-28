from __future__ import annotations

import json
import re
import unicodedata
from pathlib import Path

from .model import digest, slug, text


def key(value: str) -> str:
    normalized = unicodedata.normalize("NFKC", text(value)).upper()
    # NFKC + upper turns Greek alpha/beta into their uppercase code points;
    # handle both cases before stripping non-ASCII characters so aliases such
    # as ``AChR α`` remain distinct from the receptor-complex alias ``AChR``.
    normalized = (
        normalized.replace("Α", "ALPHA")
        .replace("α", "ALPHA")
        .replace("Β", "BETA")
        .replace("β", "BETA")
        .replace("Γ", "GAMMA")
        .replace("γ", "GAMMA")
        .replace("Δ", "DELTA")
        .replace("δ", "DELTA")
        .replace("Ε", "EPSILON")
        .replace("ε", "EPSILON")
    )
    return re.sub(r"[^A-Z0-9]+", "", normalized)


IDENTIFIER_NAMESPACES = {
    "UNIPROT": "UNIPROT",
    "HGNC": "HGNC",
    "NCBI": "NCBI_GENE",
    "NCBI_GENE": "NCBI_GENE",
    "GENEID": "NCBI_GENE",
    "ENTREZ": "NCBI_GENE",
}


def _organism_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", text(value).casefold())


def _organism_matches(entity: dict, organism: str) -> bool:
    """Compare source context with a curated entity without broad taxonomic guessing."""
    source = text(organism)
    entity_organism = entity.get("organism", {}) if isinstance(entity, dict) else {}
    if not source or not entity_organism:
        return True
    source_key = _organism_key(source)
    entity_name = text(entity_organism.get("name"))
    entity_taxonomy = text(entity_organism.get("ncbi_taxonomy_id"))
    if source_key == entity_taxonomy or source_key == _organism_key(entity_name):
        return True
    aliases = {
        "human": "homosapiens",
        "humans": "homosapiens",
        "mouse": "musmusculus",
        "murine": "musmusculus",
        "rat": "rattusnorvegicus",
        "bovine": "bostaurus",
        "cow": "bostaurus",
    }
    return aliases.get(source_key, source_key) == aliases.get(
        _organism_key(entity_name), _organism_key(entity_name)
    )


def identifier_key(value: str, namespace: str = "") -> str:
    """Normalize stable target identifiers while retaining their namespace."""
    raw = text(value).strip().upper()
    if not raw:
        return ""
    if "://" in raw:
        if "UNIPROT" in raw:
            match = re.search(r"(?:UNIPROTKB/|UNIPROT[:/])([A-Z0-9_-]+)", raw)
            if match:
                return f"UNIPROT:{match.group(1).rstrip('/')}"
        if "HGNC" in raw:
            match = re.search(r"HGNC[:/]?(\d+)", raw)
            if match:
                return f"HGNC:{match.group(1)}"
        if "ENTREZ" in raw or "NCBI" in raw or "/GENE/" in raw:
            match = re.search(r"(?:GENE|ENTREZ|ID)[/:]?(\d+)", raw)
            if match:
                return f"NCBI_GENE:{match.group(1)}"
    prefix, separator, identifier = raw.partition(":")
    selected = IDENTIFIER_NAMESPACES.get(prefix.replace(" ", "_"), "") if separator else ""
    if not selected:
        selected = IDENTIFIER_NAMESPACES.get(text(namespace).upper(), text(namespace).upper())
        identifier = raw
    return f"{selected}:{identifier.strip()}" if selected else identifier.strip()


class TargetResolver:
    """Resolve stable source identifiers first, then curated synonym overrides."""

    def __init__(self, aliases_path: Path, entities_path: Path | None = None):
        groups = json.loads(aliases_path.read_text(encoding="utf-8"))
        self.lookup: dict[str, str] = {}
        self.aliases: dict[str, set[str]] = {}
        self.identifier_lookup: dict[str, str] = {}
        self.identifier_bare_lookup: dict[str, str | None] = {}
        self.identifier_values: dict[str, set[str]] = {}
        for canonical, aliases in groups.items():
            self.aliases.setdefault(canonical, set()).update([canonical, *aliases])
            for alias in [canonical, *aliases]:
                normalized = key(alias)
                previous = self.lookup.get(normalized)
                if previous and previous != canonical:
                    raise ValueError(
                        "target alias collision: "
                        f"{alias!r} normalizes to {normalized!r} for both "
                        f"{previous!r} and {canonical!r}"
                    )
                self.lookup[normalized] = canonical
        entity_path = entities_path or aliases_path.with_name("target_entities.json")
        self.entities = (
            json.loads(entity_path.read_text(encoding="utf-8")) if entity_path.exists() else {}
        )
        for canonical, entity in self.entities.items():
            for namespace, values in (entity.get("identifiers", {}) or {}).items():
                values = values if isinstance(values, list) else [values]
                for value in values:
                    qualified = identifier_key(str(value), namespace)
                    if not qualified:
                        continue
                    previous = self.identifier_lookup.get(qualified)
                    if previous and previous != canonical:
                        raise ValueError(
                            f"target identifier collision: {qualified!r} for "
                            f"{previous!r} and {canonical!r}"
                        )
                    self.identifier_lookup[qualified] = canonical
                    self.identifier_values.setdefault(canonical, set()).add(qualified)
                    bare = qualified.split(":", 1)[-1]
                    previous_bare = self.identifier_bare_lookup.get(bare)
                    if previous_bare is None and bare in self.identifier_bare_lookup:
                        continue
                    if previous_bare and previous_bare != canonical:
                        self.identifier_bare_lookup[bare] = None
                    else:
                        self.identifier_bare_lookup[bare] = canonical

    def entity(self, canonical: str) -> dict:
        entity = dict(self.entities.get(canonical, {}))
        if entity:
            entity["mapping_provenance"] = {
                "source": "PAIRS manually curated target entity map",
                "scope": "exact canonical target only",
                "verified_on": "2026-08-18",
            }
        return entity

    def resolve(
        self, raw: str, external_id: str = "", organism: str = ""
    ) -> tuple[str, str, list[str]]:
        raw = text(raw)
        external_id = text(external_id)
        if not raw and not external_id:
            return "", "", []
        qualified = identifier_key(external_id)
        canonical = self.identifier_lookup.get(qualified, "") if qualified else ""
        resolved_by_identifier = bool(canonical)
        if not canonical and raw:
            canonical = self.identifier_lookup.get(identifier_key(raw), "")
            resolved_by_identifier = bool(canonical)
        if not canonical and raw:
            canonical = self.identifier_bare_lookup.get(raw.upper(), "") or ""
            resolved_by_identifier = bool(canonical)
        if not canonical:
            canonical = self.lookup.get(key(raw), raw.strip())
        curated = canonical in self.aliases or canonical in self.entities
        aliases = sorted(self.aliases.get(canonical, {raw, canonical}), key=str.casefold)
        for value in [raw, external_id]:
            if value and value not in aliases:
                aliases.append(value)
        if qualified and canonical in self.identifier_values:
            aliases.extend(
                value
                for value in sorted(self.identifier_values[canonical])
                if value not in aliases
            )
        if qualified and canonical not in self.identifier_values:
            # Preserve a stable unresolved identifier as the target identity;
            # free text remains its display label and searchable alias.
            target_id = "target:identifier-" + slug(qualified)
            return target_id, canonical or qualified, aliases
        if not resolved_by_identifier and organism and canonical in self.entities:
            if not _organism_matches(self.entities[canonical], organism):
                qualified_name = f"{raw or canonical} [{organism}]"
                target_id = "target:" + slug(qualified_name) + "-" + digest(
                    qualified_name.casefold(), length=12
                )
                aliases.extend(
                    value
                    for value in (raw, canonical, organism)
                    if value and value not in aliases
                )
                return target_id, qualified_name, sorted(set(aliases), key=str.casefold)
        if raw not in aliases:
            aliases.append(raw)
        target_id = "target:" + slug(canonical)
        if not curated:
            identity = unicodedata.normalize("NFKC", canonical).strip().casefold()
            target_id += "-" + digest(identity, length=12)
        return target_id, canonical, aliases

    def synonym_group(
        self, raw: str, external_id: str = "", organism: str = ""
    ) -> tuple[str, str, list[str]]:
        """Resolve explicit slash-delimited synonym groups used by Thera-SAbDab."""
        parts = [part.strip() for part in raw.split("/") if part.strip()]
        for part in parts:
            if key(part) in self.lookup:
                target_id, name, aliases = self.resolve(part, external_id, organism)
                return target_id, name, sorted(set(aliases + parts), key=str.casefold)
        base = parts[0] if parts else raw
        target_id, name, aliases = self.resolve(base, external_id, organism)
        return target_id, name, sorted(set(aliases + parts), key=str.casefold)

    def mention_group(self, raw: str) -> list[tuple[str, str, list[str]]]:
        """Resolve literature co-mentions independently; co-occurrence never creates aliases."""
        parts = [part.strip() for part in raw.split(";") if part.strip()]
        return [self.resolve(part) for part in parts]
