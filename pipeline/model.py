from __future__ import annotations

import hashlib
import re
import unicodedata
from dataclasses import dataclass, field, asdict
from typing import Any

NA = {"", "na", "n/a", "nan", "none", "nd", "not determined", "unknown", "-"}


class SequenceNormalizationError(ValueError):
    """Raised when a source sequence contains more than formatting whitespace."""


def text(value: Any) -> str:
    if value is None:
        return ""
    s = str(value).strip()
    return "" if s.lower() in NA else s


def sequence(value: Any) -> str:
    """Normalize permitted sequence formatting without deleting data.

    Whitespace is formatting; every other unexpected character is rejected so
    malformed or multipart source values cannot silently hash as a different
    biological sequence.
    """
    raw = text(value)
    normalized = re.sub(r"\s+", "", raw).upper()
    invalid = sorted(set(re.findall(r"[^A-Z*]", normalized)))
    if invalid:
        raise SequenceNormalizationError(
            f"unexpected sequence character(s): {', '.join(repr(char) for char in invalid)}"
        )
    return normalized


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
    # Source-specific structure identity tiers.  ``structures`` remains the
    # exact-structure collection used by generic retrieval/facets; tiered
    # matches (for example Thera-SAbDab's 99% SI structures) are retained
    # separately for explicit labelling.
    structure_tiers: dict[str, list[str]] = field(default_factory=dict)
    # Construct/arm context is populated by adapters that can distinguish a
    # sequence arm from the therapeutic construct.  It is deliberately kept
    # separate from target interactions.
    construct: dict[str, Any] = field(default_factory=dict)
    arm: dict[str, Any] = field(default_factory=dict)
    therapeutic_status: str = ""
    aliases: list[str] = field(default_factory=list)
    reference: str = ""
    source_url: str = ""
    record_url: str = ""
    link_scope: str = "source_homepage"
    metadata: dict[str, Any] = field(default_factory=dict)

    def is_paired(self) -> bool:
        return bool(self.heavy and self.light)

    def is_confirmed_vhh(self) -> bool:
        """Whether this record explicitly identifies a single-domain VHH."""
        # Only explicit format/domain metadata is authoritative.  Names,
        # aliases, and organism descriptions can mention nanobodies without
        # proving that this sequence is a confirmed single-domain construct.
        descriptor = " ".join([self.format, self.metadata.get("domain_type", "")]).casefold()
        return bool(
            self.heavy
            and not self.light
            and re.search(
                r"(?:\bvhh\b|\bnb\b|nanob(?:ody|odies)|single[- ]domain|single[- ]chain\s+vhh)",
                descriptor,
            )
        )

    def identity(self) -> str:
        if self.is_paired():
            return "ab_" + digest(self.heavy, self.light)
        if self.is_confirmed_vhh():
            return "ab_" + digest("vhh", self.heavy)
        # An incomplete ordinary chain is not a complete antibody identity.
        # Scope it to the source record to prevent cross-record metadata and
        # target assignments from leaking through a shared VH/VL chain.
        chain_kind = "heavy" if self.heavy else "light" if self.light else "unsequenced"
        return "ab_" + digest("source_record", self.source, self.record_id, chain_kind, self.name)


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
    assertion_origin: str = "source"
    record_url: str = ""
    link_scope: str = "source_homepage"


def public_dict(obj: Any) -> dict[str, Any]:
    return asdict(obj)
