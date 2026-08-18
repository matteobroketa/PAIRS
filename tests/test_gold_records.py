import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _antibody(antibody_id: str) -> dict:
    shard = ROOT / "data" / "v4" / "antibodies" / f"{antibody_id[3:5]}.json"
    return json.loads(shard.read_text(encoding="utf-8"))[antibody_id]


def _names(values: list[dict]) -> set[str]:
    return {item.get("name") or item.get("target_name") for item in values}


def test_manually_checked_gold_records_preserve_biological_interpretation():
    cases = json.loads((ROOT / "tests" / "gold" / "cases.json").read_text(encoding="utf-8"))
    assert len(cases) >= 4
    assert {case["source"] for case in cases} >= {
        "therasabdab",
        "covabdab",
        "poxabdab",
        "plabdab",
    }
    for case in cases:
        assert case["checked_on"] and case["record_url"].startswith("https://")
        antibody = _antibody(case["antibody_id"])
        expected = case["expected"]
        records = {
            (record.get("source"), record.get("record_id"))
            for record in antibody.get("source_records", [])
        }
        assert (case["source"], case["source_record_id"]) in records
        chain_state = (
            "paired"
            if antibody.get("heavy") and antibody.get("light")
            else "heavy_only"
            if antibody.get("heavy")
            else "light_only"
        )
        assert chain_state == expected["chains"]
        direct = _names(antibody["direct_targets"])
        assert set(expected.get("direct_targets_include", [])) <= direct
        direct_text = " ".join(direct)
        assert all(token in direct_text for token in expected.get("direct_target_contains", []))
        assert set(expected["exact_structures_include"]) <= set(
            antibody.get("exact_structures", [])
        )
        assert set(expected["negative_targets_include"]) <= _names(
            antibody.get("negative_evidence", [])
        )
        assert set(expected["literature_targets_include"]) <= _names(
            antibody.get("literature_mentions", [])
        )
