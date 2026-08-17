from __future__ import annotations

import hashlib
import re
import unicodedata
from dataclasses import dataclass, field, asdict
from typing import Any

NA = {"", "na", "n/a", "nan", "none", "nd", "not determined", "unknown", "-"}


def text(value: Any) -> str:
    if value is None:
        return ""
    s = str(value).strip()
    return "" if s.lower() in NA else s


def sequence(value: Any) -> str:
    s = text(value).upper()
    return re.sub(r"[^A-Z*]", "", s)


def split_values(value: Any, separators: str = r"[;|]") -> list[str]:
    s = text(value)
    if not s:
        return []
    return [x.strip() for x in re.split(separators, s) if x.strip() and x.strip().lower() not in NA]


def slug(value: str) -> str:
    value = unicodedata.normalize("NFKD", value).encode("ascii", "ignore").decode("ascii")
    value = value.lower().replace("α", "alpha")
    value = re.sub(r"[^a-z0-9]+", "-", value).strip("-")
    return value[:120] or "unknown"


def digest(*parts: str, length: int = 20) -> str:
    raw = "\x1f".join(parts).encode("utf-8", "replace")
    return hashlib.sha256(raw).hexdigest()[:length]


@dataclass
class AntibodyObservation:
    source: str
    record_id: str
    name: str
    heavy: str = ""
    light: str = ""
    cdrh3: str = ""
    cdrl3: str = ""
    organism: str = ""
    format: str = ""
    heavy_v: str = ""
    heavy_j: str = ""
    light_v: str = ""
    light_j: str = ""
    structures: list[str] = field(default_factory=list)
    therapeutic_status: str = ""
    aliases: list[str] = field(default_factory=list)
    reference: str = ""
    source_url: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    def identity(self) -> str:
        if self.heavy or self.light:
            return "ab_" + digest(self.heavy, self.light)
        return "ab_" + digest(self.source, self.record_id, self.name)


@dataclass
class InteractionObservation:
    antibody_id: str
    source: str
    source_record_id: str
    target_raw: str
    relationship: str
    evidence: str
    reference: str = ""
    epitope: str = ""
    assay: str = ""
    note: str = ""


def public_dict(obj: Any) -> dict[str, Any]:
    return asdict(obj)
