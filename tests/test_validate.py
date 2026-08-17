import csv
import json
from pathlib import Path

from pipeline.build import SCHEMA_VERSION, compile_data, load_sources
from pipeline.validate import validate


def _tiny_dataset(tmp_path: Path) -> Path:
    source = tmp_path / "thera.csv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Therapeutic",
                "HeavySequence",
                "LightSequence",
                "Target",
                "Format",
                "Est. Status",
                "Highest_Clin_Trial (Feb '25)",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "Therapeutic": "Integritymab",
                "HeavySequence": "EVQLVESGGGLVQPGGSLRLS",
                "LightSequence": "DIQMTQSPSSLSASVGDRVT",
                "Target": "ERBB2/HER2",
                "Format": "Whole mAb",
                "Est. Status": "Active",
                "Highest_Clin_Trial (Feb '25)": "Phase-III",
            }
        )

    out = tmp_path / "out"
    stats = compile_data({"therasabdab": source}, load_sources(), out)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "app_version": "test",
        "stats": stats,
        "sources_expected": 1,
        "sources_ok": 1,
    }
    (out / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return out


def test_validate_accepts_complete_generated_dataset(tmp_path: Path):
    out = _tiny_dataset(tmp_path)
    assert validate(out) == []


def test_validate_rejects_target_reference_to_missing_antibody(tmp_path: Path):
    out = _tiny_dataset(tmp_path)
    target = json.loads((out / "targets.json").read_text(encoding="utf-8"))[0]
    page = out / "targets" / target["dir"] / "page-001.json"
    payload = json.loads(page.read_text(encoding="utf-8"))
    payload[0]["antibody"]["id"] = "ab_missing"
    page.write_text(json.dumps(payload), encoding="utf-8")

    errors = validate(out)
    assert any("references missing antibody ab_missing" in error for error in errors)


def test_validate_rejects_non_direct_relationship_on_primary_page(tmp_path: Path):
    out = _tiny_dataset(tmp_path)
    target = json.loads((out / "targets.json").read_text(encoding="utf-8"))[0]
    page = out / "targets" / target["dir"] / "page-001.json"
    payload = json.loads(page.read_text(encoding="utf-8"))
    payload[0]["interactions"][0]["relationship"] = "mentioned_with"
    page.write_text(json.dumps(payload), encoding="utf-8")

    errors = validate(out)
    assert any("unexpected relationship in pages" in error for error in errors)
