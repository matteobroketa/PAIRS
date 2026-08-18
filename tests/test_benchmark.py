import json

from pipeline.benchmark import benchmark, controlled_variant
from pipeline.build import write_similarity_indexes


def test_controlled_variant_has_deterministic_identity():
    sequence = "ACDEFGHIKLMNPQRSTVWY" * 5
    variant = controlled_variant(sequence, 90)
    assert len(variant) == len(sequence)
    assert sum(left == right for left, right in zip(sequence, variant, strict=True)) == 90


def test_similarity_benchmark_measures_candidate_recall(tmp_path):
    sequence = "ACDEFGHIKLMNPQRSTVWY" * 5
    antibodies = {
        "ab_00reference0000000000": {"heavy": sequence, "light": sequence[::-1]},
        "ab_01decoy0000000000000": {"heavy": sequence[20:] + sequence[:20], "light": sequence},
    }
    (tmp_path / "antibodies").mkdir()
    (tmp_path / "antibodies" / "00.json").write_text(
        json.dumps({"ab_00reference0000000000": antibodies["ab_00reference0000000000"]}),
        encoding="utf-8",
    )
    (tmp_path / "antibodies" / "01.json").write_text(
        json.dumps({"ab_01decoy0000000000000": antibodies["ab_01decoy0000000000000"]}),
        encoding="utf-8",
    )
    write_similarity_indexes(antibodies, tmp_path)
    result = benchmark(tmp_path, samples=2)
    assert result["chains"]["heavy"]["levels"]["100"]["recall"] == 1.0
    assert result["chains"]["light"]["levels"]["100"]["recall"] == 1.0
