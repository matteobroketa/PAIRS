from pipeline.sources import _cov_target, _pox_entity


def test_cov_parent_target_is_labelled_as_derived_hierarchy():
    targets = dict(_cov_target("SARS-CoV-2", "RBD"))
    assert targets["SARS-CoV-2 RBD"] == "source_epitope"
    assert targets["SARS-CoV-2 Spike"] == "derived_hierarchy"


def test_pox_target_retains_strain_location_and_protein():
    target, raw = _pox_entity("VACV_WR_Copenhagen_A27")
    assert target == "Vaccinia virus WR Copenhagen A27"
    assert raw == "VACV_WR_Copenhagen_A27"
