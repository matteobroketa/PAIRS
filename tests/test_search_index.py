import json
from pathlib import Path

from pipeline.search import build_exact_antibody_index, normalize_search_term, search_bucket


def test_exact_index_keeps_alias_after_compact_suggestion_boundary():
    antibody = {
        "id": "ab_a1234567890123456789",
        "name": "Examplemab",
        "aliases": [*(f"alias{index:02d}" for index in range(1, 20)), "trastuzumab"],
    }

    exact = build_exact_antibody_index({antibody["id"]: antibody})

    normalized = normalize_search_term("trastuzumab")
    assert exact[search_bucket(normalized)][normalized] == [antibody["id"]]
    assert antibody["aliases"][:12][-1] != "trastuzumab"


def test_exact_index_preserves_duplicate_name_entities():
    antibodies = {
        "ab_a1234567890123456789": {"name": "ExampleAb", "aliases": []},
        "ab_b1234567890123456789": {"name": "ExampleAb", "aliases": []},
    }

    exact = build_exact_antibody_index(antibodies)

    normalized = normalize_search_term("ExampleAb")
    assert exact[search_bucket(normalized)][normalized] == sorted(antibodies)


def test_generated_exact_index_is_a_term_to_id_map(tmp_path: Path):
    from pipeline.build import write_indexes

    antibodies = {
        "ab_a1234567890123456789": {
            "id": "ab_a1234567890123456789",
            "name": "Roundtripmab",
            "aliases": ["Round Trip Alias"],
            "sources": [],
            "structures": [],
        }
    }
    write_indexes(antibodies, [], tmp_path / "out", {})

    payload = json.loads((tmp_path / "out" / "antibody-exact" / "ro.json").read_text())
    assert payload["roundtripmab"] == ["ab_a1234567890123456789"]
    assert payload["round trip alias"] == ["ab_a1234567890123456789"]
