"""Shared normalization and indexing rules for textual PAIRS search."""

from __future__ import annotations

import unicodedata
from collections import defaultdict
from collections.abc import Iterable


def normalize_search_term(value: object) -> str:
    """Return the browser-compatible normalized textual search term.

    PAIRS deliberately treats punctuation runs as separators so established
    forms such as ``PD-L1`` and ``PDL1`` remain equivalent.  This is limited
    to textual names and aliases; it is not used for biological sequences.
    """

    value = "" if value is None else str(value)
    value = unicodedata.normalize("NFKD", value).lower()
    value = "".join(
        character if character.isascii() and character.isalnum() else " "
        for character in value
    )
    return " ".join(value.split())


def search_bucket(normalized_term: str) -> str:
    """Return the two-character lazy-load bucket for a normalized term."""

    compact = normalized_term.replace(" ", "")
    if not compact:
        return "__"
    return compact[:2] if len(compact) >= 2 else f"{compact}_"


def antibody_search_terms(antibody: dict) -> list[str]:
    """Return every retained source-derived term supported by antibody search."""

    values: list[str] = [antibody.get("name", ""), *(antibody.get("aliases", []) or [])]
    # These optional fields are not currently emitted by the public adapters,
    # but keep the exact index contract lossless if a source exposes them.
    for field in ("therapeutic_name", "therapeutic_names", "public_identifiers"):
        value = antibody.get(field, "")
        values.extend(value if isinstance(value, list) else [value])
    return [value for value in values if isinstance(value, str) and value.strip()]


def build_exact_antibody_index(antibodies: dict[str, dict]) -> dict[str, dict[str, list[str]]]:
    """Build ``bucket -> normalized term -> antibody IDs`` without truncation."""

    buckets: dict[str, defaultdict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for antibody_id, antibody in antibodies.items():
        for raw_term in antibody_search_terms(antibody):
            normalized = normalize_search_term(raw_term)
            if normalized:
                buckets[search_bucket(normalized)][normalized].add(antibody_id)
    return {
        bucket: {
            term: sorted(ids)
            for term, ids in sorted(terms.items())
        }
        for bucket, terms in sorted(buckets.items())
    }


def iter_exact_antibody_terms(antibodies: dict[str, dict]) -> Iterable[tuple[str, str, str]]:
    """Yield ``(antibody_id, raw_term, normalized_term)`` for round-trip checks."""

    for antibody_id, antibody in antibodies.items():
        for raw_term in antibody_search_terms(antibody):
            normalized = normalize_search_term(raw_term)
            if normalized:
                yield antibody_id, raw_term, normalized
