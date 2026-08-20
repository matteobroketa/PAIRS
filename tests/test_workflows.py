from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_pages_deploy_never_rebuilds_production_data():
    workflow = (ROOT / ".github" / "workflows" / "deploy-pages.yml").read_text(
        encoding="utf-8"
    )
    assert "pipeline.build" not in workflow
    assert "actions/cache" not in workflow
    assert "migration-cache" not in workflow
    assert "pairs-production-v4" in workflow
    assert "refusing capped starter data" in workflow


def test_full_build_triggers_for_pipeline_and_config_changes():
    workflow = (ROOT / ".github" / "workflows" / "build-data.yml").read_text(
        encoding="utf-8"
    )
    assert "push:" in workflow
    assert '"pipeline/**"' in workflow
    assert '"config/**"' in workflow
    assert "schedule:" in workflow
    assert "workflow_dispatch:" in workflow
    assert "pipeline.build --output production-data/v4" in workflow
    assert "pairs-production-v4" in workflow
