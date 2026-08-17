from pipeline.model import AntibodyObservation, sequence


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
