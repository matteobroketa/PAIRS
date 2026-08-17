import pytest

from pipeline.model import AntibodyObservation, SequenceNormalizationError, sequence


def test_sequence_normalization():
    assert sequence(" EVQL VES\nGG ") == "EVQLVESGG"


def test_exact_pair_identity_is_source_independent():
    a=AntibodyObservation(source="a",record_id="1",name="x",heavy="AAA",light="BBB")
    b=AntibodyObservation(source="b",record_id="2",name="y",heavy="AAA",light="BBB")
    assert a.identity() == b.identity()


def test_unsequenced_records_do_not_cross_merge():
    a=AntibodyObservation(source="a",record_id="1",name="same")
    b=AntibodyObservation(source="b",record_id="2",name="same")
    assert a.identity() != b.identity()


def test_incomplete_heavy_and_light_records_are_source_record_scoped():
    heavy_a = AntibodyObservation(source="covabdab", record_id="1", name="a", heavy="AAA")
    heavy_b = AntibodyObservation(source="covabdab", record_id="2", name="b", heavy="AAA")
    light_a = AntibodyObservation(source="covabdab", record_id="3", name="c", light="BBB")
    light_b = AntibodyObservation(source="covabdab", record_id="4", name="d", light="BBB")
    assert heavy_a.identity() != heavy_b.identity()
    assert light_a.identity() != light_b.identity()


def test_confirmed_vhh_can_deduplicate_by_single_domain():
    a = AntibodyObservation(source="a", record_id="1", name="x", heavy="AAA", format="VHH")
    b = AntibodyObservation(source="b", record_id="2", name="y", heavy="AAA", format="nanobody")
    assert a.identity() == b.identity()


def test_sequence_rejects_unexpected_characters_instead_of_deleting_them():
    with pytest.raises(SequenceNormalizationError):
        sequence("EVQL- VESGG")
