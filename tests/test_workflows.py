from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_pages_deploy_never_rebuilds_production_data():
    workflow = (ROOT / ".github" / "workflows" / "deploy-pages.yml").read_text(
        encoding="utf-8"
    )
    assert "pipeline.build" not in workflow
    assert "pairs-production-v4" in workflow
    assert "refusing capped starter data" in workflow


def test_full_build_is_isolated_to_scheduled_or_manual_workflow():
    workflow = (ROOT / ".github" / "workflows" / "build-data.yml").read_text(
        encoding="utf-8"
    )
    assert "push:" not in workflow
    assert "schedule:" in workflow
    assert "workflow_dispatch:" in workflow
    assert "pipeline.build --output production-data/v4" in workflow
    assert "pairs-production-v4" in workflow
