from __future__ import annotations

import json
import re
import unicodedata
from pathlib import Path

from .model import digest, slug, text


def key(value: str) -> str:
    normalized = unicodedata.normalize("NFKC", text(value)).upper()
    normalized = normalized.replace("α", "ALPHA").replace("Β", "BETA").replace("β", "BETA")
    return re.sub(r"[^A-Z0-9]+", "", normalized)


class TargetResolver:
    """Conservative target resolver backed only by explicitly curated synonym groups."""

    def __init__(self, aliases_path: Path, entities_path: Path | None = None):
        groups = json.loads(aliases_path.read_text(encoding="utf-8"))
        self.lookup: dict[str, str] = {}
        self.aliases: dict[str, set[str]] = {}
        for canonical, aliases in groups.items():
            self.aliases.setdefault(canonical, set()).update([canonical, *aliases])
            for alias in [canonical, *aliases]:
                self.lookup[key(alias)] = canonical
        entity_path = entities_path or aliases_path.with_name("target_entities.json")
        self.entities = (
            json.loads(entity_path.read_text(encoding="utf-8")) if entity_path.exists() else {}
        )

    def entity(self, canonical: str) -> dict:
        entity = dict(self.entities.get(canonical, {}))
        if entity:
            entity["mapping_provenance"] = {
                "source": "PAIRS manually curated target entity map",
                "scope": "exact canonical target only",
                "verified_on": "2026-08-18",
            }
        return entity

    def resolve(self, raw: str) -> tuple[str, str, list[str]]:
        raw = text(raw)
        if not raw:
            return "", "", []
        curated = key(raw) in self.lookup
        canonical = self.lookup.get(key(raw), raw.strip())
        aliases = sorted(self.aliases.get(canonical, {raw, canonical}), key=str.casefold)
        if raw not in aliases:
            aliases.append(raw)
        target_id = "target:" + slug(canonical)
        if not curated:
            identity = unicodedata.normalize("NFKC", canonical).strip().casefold()
            target_id += "-" + digest(identity, length=12)
        return target_id, canonical, aliases

    def synonym_group(self, raw: str) -> tuple[str, str, list[str]]:
        """Resolve explicit slash-delimited synonym groups used by Thera-SAbDab."""
        parts = [part.strip() for part in raw.split("/") if part.strip()]
        for part in parts:
            if key(part) in self.lookup:
                target_id, name, aliases = self.resolve(part)
                return target_id, name, sorted(set(aliases + parts), key=str.casefold)
        base = parts[0] if parts else raw
        target_id, name, aliases = self.resolve(base)
        return target_id, name, sorted(set(aliases + parts), key=str.casefold)

    def mention_group(self, raw: str) -> list[tuple[str, str, list[str]]]:
        """Resolve literature co-mentions independently; co-occurrence never creates aliases."""
        parts = [part.strip() for part in raw.split(";") if part.strip()]
        return [self.resolve(part) for part in parts]
