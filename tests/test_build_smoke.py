import csv
import json
from pathlib import Path

from pipeline.build import compile_data, load_sources


def test_tiny_thera_build(tmp_path: Path):
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
                "Therapeutic": "Examplemab",
                "HeavySequence": "AAA",
                "LightSequence": "BBB",
                "Target": "ERBB2/CD340/HER2",
                "Format": "Whole mAb",
                "Est. Status": "Active",
                "Highest_Clin_Trial (Feb '25)": "Phase-III",
            }
        )
    config = load_sources()
    stats = compile_data({"therasabdab": source}, config, tmp_path / "out")
    assert stats["antibodies"] == 1
    targets = json.loads((tmp_path / "out" / "targets.json").read_text())
    assert any(target["name"] == "ERBB2" for target in targets)
    assert (tmp_path / "out" / "sequence-search").exists()


def test_plabdab_comentions_fan_out_without_alias_merging(tmp_path: Path):
    source = tmp_path / "plabdab.csv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["ID", "heavy_sequence", "light_sequence", "targets_mentioned"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "ID": "paper-antibody",
                "heavy_sequence": "EVQLVESGG",
                "light_sequence": "DIQMTQSPSS",
                "targets_mentioned": "ACP; ACP1; NDUFAB1; SDAP; FASN2A",
            }
        )
    config = load_sources()
    stats = compile_data({"plabdab": source}, config, tmp_path / "out")
    assert stats["targets"] == 5
    assert stats["interactions"] == 5
    targets = json.loads((tmp_path / "out" / "targets.json").read_text())
    by_name = {target["name"]: target for target in targets}
    assert set(by_name) == {"ACP", "ACP1", "NDUFAB1", "SDAP", "FASN2A"}
    assert "ACP1" not in by_name["ACP"]["aliases"]
    assert "ACP" not in by_name["ACP1"]["aliases"]
