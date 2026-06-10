from click.testing import CliRunner
import pytest

from metatrawl import cli
from metatrawl import healthcheck


def test_healthcheck_reports_all_requested_dependencies():
    def executable_probe(command: str) -> tuple[bool, str]:
        return True, f"/usr/bin/{command}"

    def package_probe(import_name: str, package_name: str | None) -> tuple[bool, str]:
        return True, f"installed ({import_name})"

    checks = healthcheck.collect_checks(
        executable_probe=executable_probe,
        package_probe=package_probe,
    )

    assert [check.name for check in checks] == [
        "zipstrain",
        "sylph",
        "samtools",
        "bowtie2",
        "prefetch",
        "fasterq-dump",
        "datasets",
        "prodigal",
        "torch",
        "h5py",
    ]
    assert {check.status for check in checks} == {"ok"}


def test_metatrawl_test_renders_dependency_report(monkeypatch):
    monkeypatch.setattr(
        healthcheck,
        "collect_checks",
        lambda: [
            healthcheck.Check("zipstrain", "ok", "/usr/bin/zipstrain"),
            healthcheck.Check("sylph", "missing", "not found in PATH"),
        ],
    )

    result = CliRunner().invoke(cli.cli, ["test"])

    assert result.exit_code == 0
    assert "MetaTrawl health check" in result.output
    assert "zipstrain" in result.output
    assert "sylph" in result.output
    assert "MISSING" in result.output


def test_assert_dependencies_fails_with_clear_missing_list():
    checks = [
        healthcheck.Check("zipstrain", "ok", "/usr/bin/zipstrain"),
        healthcheck.Check("sylph", "missing", "not found in PATH"),
    ]

    with pytest.raises(healthcheck.DependencyCheckError) as exc_info:
        healthcheck.assert_dependencies(checks)

    assert "sylph" in str(exc_info.value)
    assert "not found in PATH" in str(exc_info.value)


def test_metatrawl_check_fails_when_dependency_missing(monkeypatch):
    monkeypatch.setattr(
        healthcheck,
        "collect_checks",
        lambda: [
            healthcheck.Check("zipstrain", "ok", "/usr/bin/zipstrain"),
            healthcheck.Check("sylph", "missing", "not found in PATH"),
        ],
    )

    result = CliRunner().invoke(cli.cli, ["check"])

    assert result.exit_code != 0
    assert "MetaTrawl health check" in result.output
    assert "Missing required MetaTrawl dependencies" in result.output
