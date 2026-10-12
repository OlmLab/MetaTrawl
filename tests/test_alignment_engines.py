"""Alignment protocol, configuration, provenance, and worker handoff tests."""

from dataclasses import replace
import hashlib
import json
from pathlib import Path
import random
import shlex
import shutil
import subprocess

from click.testing import CliRunner
import pytest

from metatrawl import alignment, cache, cli, db, provenance, workflows
from metatrawl.config import ProfileConfig, WorkflowConfig, load_workflow_config
from metatrawl.execution import WorkflowRuntime
from metatrawl.logging import WorkflowLogger
from test_provenance import bundle


def test_alignment_engine_is_abstract():
    with pytest.raises(TypeError, match="abstract"):
        alignment.AlignmentEngine()


@pytest.mark.parametrize("name", ["bowtie2", "bwa-mem3"])
@pytest.mark.parametrize("paired", [False, True])
def test_engine_commands_and_threads(name, paired):
    engine = alignment.get_alignment_engine(name)
    ref = Path("reference.fna")
    reads = [Path("R1.fastq.gz")]
    if paired:
        reads.append(Path("R2.fastq.gz"))
    if name == "bowtie2":
        assert engine.index_command(ref, 7) == ["bowtie2-build", "--threads", "7", str(ref), str(ref)]
        inputs = ["-1", "R1.fastq.gz", "-2", "R2.fastq.gz"] if paired else ["-U", "R1.fastq.gz"]
        expected = ["bowtie2", "-x", str(ref), *inputs, "--threads", "9"]
    else:
        assert engine.index_command(ref, 7) == ["bwa-mem3", "index", "-t", "7", str(ref)]
        expected = ["bwa-mem3", "mem", "-t", "9", str(ref), *map(str, reads)]
    assert engine.alignment_command(ref, reads, 9) == expected
    assert "--fast" not in expected
    command = engine.bam_command(ref, reads, 9, Path("sample.bam.tmp"))
    argv = shlex.split(command)
    assert argv[:4] == ["bash", "-o", "pipefail", "-c"]
    assert shlex.join(expected) in argv[4]
    assert "samtools view -bS -F 4 -" in argv[4]
    assert "samtools sort -@ 9 -o sample.bam.tmp -" in argv[4]


@pytest.mark.parametrize("name", ["bowtie2", "bwa-mem3"])
def test_alignment_quotes_paths(name):
    engine = alignment.get_alignment_engine(name)
    ref = Path("odd path's/reference.fna")
    reads = [Path("odd path's/read 1.fq.gz"), Path("odd path's/read 2.fq.gz")]
    bam = Path("odd path's/result.bam.tmp")
    script = shlex.split(engine.bam_command(ref, reads, 2, bam))[4]
    assert shlex.join(engine.alignment_command(ref, reads, 2)) in script
    assert f"-o {shlex.quote(str(bam))} -" in script


@pytest.mark.parametrize("large", [False, True])
def test_bowtie_index_small_large_and_partial(tmp_path, large):
    engine = alignment.get_alignment_engine()
    ref = tmp_path / "ref.fna"
    paths = engine.index_paths(ref)[6:] if large else engine.index_paths(ref)[:6]
    assert not engine.index_complete(ref)
    for path in paths[:-1]:
        path.write_text("index")
    assert not engine.index_complete(ref)
    paths[-1].touch()
    assert not engine.index_complete(ref)
    paths[-1].write_text("index")
    assert engine.index_complete(ref)


def test_bwa_index_and_optional_unpacked_cleanup(tmp_path):
    engine = alignment.get_alignment_engine("bwa-mem3")
    ref = tmp_path / "ref.fna"
    assert [path.name for path in engine.index_paths(ref)] == [
        "ref.fna.amb", "ref.fna.ann", "ref.fna.pac", "ref.fna.bwt.2bit.64",
    ]
    for path in engine.index_paths(ref)[:-1]:
        path.write_text("index")
    assert not engine.index_complete(ref)
    engine.index_paths(ref)[-1].write_text("index")
    assert engine.index_complete(ref)
    assert Path(f"{ref}.0123") in engine.cleanup_paths(ref)


@pytest.mark.parametrize("suffix", ["toml", "json"])
@pytest.mark.parametrize("name", ["bowtie2", "bwa-mem3"])
def test_engine_configuration(tmp_path, suffix, name):
    path = tmp_path / f"workflow.{suffix}"
    path.write_text(
        f'[profile]\nalignment_engine = "{name}"\n' if suffix == "toml"
        else json.dumps({"profile": {"alignment_engine": name}})
    )
    config = load_workflow_config(path, threads=2, sample_count=1)
    assert config.profile.alignment_engine == name
    assert config.profile.backend == "python"


@pytest.mark.parametrize("name", [None, "", "   ", 1, False, "no-such-aligner"])
def test_invalid_engine_configuration(tmp_path, name):
    path = tmp_path / "workflow.json"
    path.write_text(json.dumps({"profile": {"alignment_engine": name}}))
    with pytest.raises(ValueError, match="engine"):
        load_workflow_config(path, threads=2, sample_count=1)


@pytest.mark.parametrize("command", ["sync-profile", "profile-sra"])
@pytest.mark.parametrize("configured,override,expected", [
    (None, None, "bowtie2"), ("bwa-mem3", None, "bwa-mem3"),
    ("bowtie2", "bwa-mem3", "bwa-mem3"), ("bwa-mem3", "bowtie2", "bowtie2"),
])
def test_cli_engine_precedence(tmp_path, monkeypatch, command, configured, override, expected):
    configs = []

    def sync(**kwargs):
        configs.append(kwargs["workflow_config"])
        return workflows.SyncSummary(0, 0, 0, 0, 0)

    def profile(**kwargs):
        configs.append(kwargs["workflow_config"])
        return {}

    monkeypatch.setattr(workflows, "sync_remaining_profiles", sync)
    monkeypatch.setattr(workflows, "profile_sra_runs", profile)
    args = [command, "--db", str(tmp_path / "project.duckdb"), "--cache-dir", "cache",
            "--scratch-dir", "scratch", "--output-dir", "outputs", "--backend", "rust_profiler"]
    if command == "profile-sra":
        remaining = tmp_path / "remaining.csv"
        remaining.write_text("run_id\nSRR1\n")
        args.extend(["--remaining-csv", str(remaining)])
    if configured is not None:
        path = tmp_path / "workflow.toml"
        path.write_text(f'[profile]\nalignment_engine = "{configured}"\n')
        args.extend(["--workflow-config", str(path)])
    if override is not None:
        args.extend(["--alignment-engine", override])
    result = CliRunner().invoke(cli.cli, args)
    assert result.exit_code == 0, result.output
    assert configs[0].profile.alignment_engine == expected
    assert configs[0].profile.backend == "rust_profiler"


def test_default_contract_unchanged_and_bwa_import_requires_acknowledgement(tmp_path):
    default = provenance.profiling_contract(ProfileConfig())
    assert set(default["profile"]) == {"min_mapq", "min_baseq", "min_freq", "min_read_ani", "read_inclusion"}
    assert default["alignment"] == {
        "aligner": "bowtie2", "preset": "default-sensitive-end-to-end", "samtools_exclude_flags": 4,
    }
    bwa = replace(ProfileConfig(), alignment_engine="bwa-mem3")
    contract = provenance.profiling_contract(bwa)
    assert contract["profile"] == default["profile"]
    assert "alignment.aligner" in "; ".join(provenance.differences(default, contract))
    with db.connect(tmp_path / "project.duckdb") as conn:
        db.add_runs(conn, ["S1", "S2"])
        db.import_profile_bundle(conn, bundle(tmp_path, "S1"))
        second = bundle(tmp_path, "S2", config=bwa)
        with pytest.raises(ValueError, match="alignment.aligner"):
            db.import_profile_bundle(conn, second)
        assert conn.execute("SELECT count(*) FROM profiles").fetchone()[0] == 1
        db.import_profile_bundle(conn, second, allow_incompatible_profiles=True)
        assert conn.execute("SELECT count(*) FROM profiling_contracts").fetchone()[0] == 2


def test_alignment_checkpoint_cannot_switch_engines(tmp_path):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    (scratch / "profiling_settings.json").write_text(provenance.canonical(provenance.profiling_contract(ProfileConfig())))
    reference = cache.PreparedReference(Path("ref.fna"), Path("genes.fna"), Path("ref.stb"))
    with pytest.raises(ValueError, match="alignment.aligner"):
        workflows._run_alignment_and_profile(
            run_id="S1", sample_scratch=scratch, reference=reference,
            output_dir=tmp_path / "outputs", logger=WorkflowLogger(),
            profile_config=replace(ProfileConfig(), alignment_engine="bwa-mem3"),
        )


def test_missing_bwa_errors_in_execution_environment(tmp_path):
    config = WorkflowConfig.legacy(threads=2, sample_count=1)
    config = replace(config, stages={**config.stages, "bowtie_build": replace(
        config.stage("bowtie_build"), environment={"PATH": "/usr/bin:/bin"})})
    runtime = WorkflowRuntime(config, state_dir=tmp_path, logger=WorkflowLogger())
    with pytest.raises(RuntimeError, match="requires 'bwa-mem3'.*Install it separately"):
        alignment.get_alignment_engine("bwa-mem3").prepare_index(runtime, tmp_path / "ref.fna", sample="S1")


@pytest.mark.parametrize("name", ["bowtie2", "bwa-mem3"])
def test_dependency_check_does_not_require_unused_bowtie(tmp_path, monkeypatch, name):
    selected = []
    monkeypatch.setattr(workflows.healthcheck, "collect_checks", lambda **kwargs: selected.append(kwargs["executables"]) or [])
    monkeypatch.setattr(workflows.healthcheck, "assert_dependencies", lambda checks: checks)
    config = replace(WorkflowConfig.legacy(threads=1, sample_count=1),
                     profile=replace(ProfileConfig(), alignment_engine=name))
    summary = workflows.sync_remaining_profiles(
        db_file=tmp_path / "project.duckdb", cache_dir=tmp_path / "cache",
        scratch_dir=tmp_path / "scratch", output_dir=tmp_path / "outputs", workflow_config=config,
    )
    assert summary.requested == 0
    assert ("bowtie2" in selected[0]) == (name == "bowtie2")
    assert "samtools" in selected[0]


def test_failed_aligner_not_hidden_by_successful_sort(tmp_path):
    tools = tmp_path / "bin"
    tools.mkdir()
    (tools / "bwa-mem3").write_text("#!/bin/sh\nexit 23\n")
    (tools / "samtools").write_text(
        '#!/bin/sh\nif [ "$1" = view ]; then cat; exit 0; fi\n'
        'while [ "$#" -gt 0 ]; do if [ "$1" = -o ]; then shift; output=$1; fi; shift; done\n'
        'cat >/dev/null\nprintf partial > "$output"\n'
    )
    for path in tools.iterdir():
        path.chmod(0o755)
    config = WorkflowConfig.legacy(threads=2, sample_count=1)
    config = replace(config, stages={**config.stages, "alignment": replace(
        config.stage("alignment"), environment={"PATH": f"{tools}:/usr/bin:/bin"})})
    runtime = WorkflowRuntime(config, state_dir=tmp_path, logger=WorkflowLogger())
    bam = tmp_path / "partial.bam"
    command = alignment.get_alignment_engine("bwa-mem3").bam_command(Path("ref.fna"), [Path("read.fq")], 2, bam)
    with pytest.raises(RuntimeError, match="command failed"):
        runtime.run_shell("alignment", command, sample="S1")
    assert bam.read_text() == "partial"


def test_fixed_reference_indexes_cached_separately_by_engine(tmp_path):
    ref, genes, stb = (tmp_path / name for name in ("ref.fna", "genes.fna", "ref.stb"))
    ref.write_text(">chr\nACGT\n")
    genes.write_text(">chr_1 # 1 # 4 # 1 # ID=1\nACGT\n")
    stb.write_text("chr\tgenome\n")
    calls = []

    class Runtime:
        def threads(self, stage):
            return 3

        def run(self, stage, command, **kwargs):
            calls.append(command)
            if command[0] == "zipstrain":
                output = Path(command[command.index("--output-dir") + 1])
                output.mkdir()
                for name in ("reference.fasta", "genomes_bed_file.bed", "gene_range_table.tsv", "null_model.parquet", "profiling_contract.json"):
                    (output / name).write_text("fixture")
            elif command[:2] == ["samtools", "faidx"]:
                Path(f"{command[2]}.fai").write_text("fixture")
            else:
                for path in alignment.get_alignment_engine().index_paths(Path(command[-1]))[:6]:
                    path.write_text("index")

        def run_shell(self, stage, command, **kwargs):
            argv = shlex.split(command.split("; exec ", 1)[1])
            calls.append(argv)
            assert argv[:4] == ["bwa-mem3", "index", "-t", "3"]
            for path in alignment.get_alignment_engine("bwa-mem3").index_paths(Path(argv[-1])):
                path.write_text("index")

    contexts = []
    for name in ("bowtie2", "bwa-mem3"):
        kwargs = dict(reference_genome=ref, reference_genome_genes=genes, reference_stb=stb,
                      cache_dir=tmp_path / "cache", runtime=Runtime(), logger=WorkflowLogger(), alignment_engine=name)
        first = workflows._prepare_fixed_reference_context(**kwargs)
        assert workflows._prepare_fixed_reference_context(**kwargs) == first
        assert alignment.get_alignment_engine(name).index_complete(first.reference.reference_fasta)
        contexts.append(first)
    assert contexts[0].reference_key != contexts[1].reference_key
    original = {"sources": {label: provenance.sha256_file(path) for label, path in {"reference": ref, "genes": genes, "stb": stb}.items()},
                "null_model": provenance.zipstrain_defaults()["null_model"]}
    assert contexts[0].reference_key == hashlib.sha256(provenance.canonical(original).encode()).hexdigest()
    assert sum(command[0] == "bowtie2-build" for command in calls) == 1
    assert sum(command[:2] == ["bwa-mem3", "index"] for command in calls) == 1


@pytest.mark.parametrize("name", ["bowtie2", "bwa-mem3"])
@pytest.mark.parametrize("execution", ["local", "slurm"])
def test_worker_engine_selection_cleanup_and_profile_handoff(tmp_path, monkeypatch, name, execution):
    sample = tmp_path / "scratch" / "S1"
    work = sample / "zipstrain_profile"
    ref_dir = sample / "reference"
    work.mkdir(parents=True)
    ref_dir.mkdir()
    ref, genes, stb = (ref_dir / path for path in ("ref.fna", "genes.fna", "ref.stb"))
    ref.write_text(">chr\nACGT\n")
    genes.write_text(">chr_1 # 1 # 4 # 1 # ID=1\nACGT\n")
    stb.write_text("chr\tgenome\n")
    for filename in ("genomes_bed_file.bed", "gene_range_table.tsv", "null_model.parquet", "reference.fasta"):
        (work / filename).write_text("fixture")
    (work / "profiling_contract.json").write_text("{}")
    (sample / "S1_1.fastq.gz").write_text("fixture")
    (sample / "S1_2.fastq.gz").write_text("fixture")
    profile = replace(ProfileConfig(), alignment_engine=name, backend="rust_profiler", min_freq=0.02)
    (sample / "profiling_settings.json").write_text(provenance.canonical(provenance.profiling_contract(profile)))
    monkeypatch.setattr(workflows, "_ensure_reference_fai", lambda **kwargs: None)
    monkeypatch.setattr(workflows, "_validate_bam", lambda *args, **kwargs: None)
    calls, scripts = [], []

    def runner(command, **kwargs):
        if isinstance(command, list) and command[0] == "sbatch":
            script = Path(command[-1]).read_text()
            scripts.append(script)
            command = script.splitlines()[-1]
        calls.append(command)
        rendered = command if isinstance(command, str) else shlex.join(command)
        if "samtools sort" in rendered:
            (sample / "S1.bam.tmp").write_text("bam")
        elif "profile-single" in rendered:
            for suffix in ("profile", "genome_stats", "gene_stats"):
                (work / f"S1_{suffix}.parquet").write_text("fixture")
        else:
            engine = alignment.get_alignment_engine(name)
            paths = engine.index_paths(ref)[:6] if name == "bowtie2" else engine.index_paths(ref)
            for path in paths:
                path.write_text("index")
        return subprocess.CompletedProcess(command, 0, stdout="123\n", stderr="")

    config = replace(WorkflowConfig.legacy(threads=6, sample_count=1), profile=profile)
    config = replace(config, stages={key: replace(value, execution=execution) for key, value in config.stages.items()})
    runtime = WorkflowRuntime(config, state_dir=tmp_path / "state", logger=WorkflowLogger(), runner=runner)
    outputs = tmp_path / "outputs"
    outputs.mkdir()
    (outputs / "S1.sylph.tsv").write_text("Genome_file\tabundance\ngenome\t1\n")
    workflows._run_alignment_and_profile(
        run_id="S1", sample_scratch=sample, reference=cache.PreparedReference(ref, genes, stb),
        output_dir=outputs, runtime=runtime, profile_config=profile, logger=WorkflowLogger(),
    )
    assert len(scripts) == (3 if execution == "slurm" else 0)
    profile_command = shlex.split(calls[-1]) if isinstance(calls[-1], str) else calls[-1]
    assert profile_command[profile_command.index("--backend") + 1] == "rust_profiler"
    assert profile_command[profile_command.index("--max-concurrency") + 1] == "6"
    assert "--alignment-engine" not in profile_command
    manifest = json.loads(provenance.manifest_path(outputs / "S1.profile.parquet").read_text())
    assert manifest["details"]["alignment_engine"] == name
    assert manifest["contract"]["alignment"]["aligner"] == name
    assert (outputs / ".S1.bundle.ready").exists()
    assert not (sample / "S1.bam").exists()
    assert not list(sample.rglob("*.fastq.gz"))
    assert not ref_dir.exists()


def test_adapter_entry_point_protocol(monkeypatch):
    class Custom(alignment.BwaMem3AlignmentEngine):
        name = "custom"

    class Entry:
        def load(self):
            return Custom

    monkeypatch.setattr(alignment, "entry_points", lambda **kwargs: [Entry()])
    assert isinstance(alignment.get_alignment_engine("custom"), alignment.AlignmentEngine)
    monkeypatch.setattr(alignment, "entry_points", lambda **kwargs: [Entry(), Entry()])
    with pytest.raises(ValueError, match="Multiple alignment"):
        alignment.get_alignment_engine("custom")
    monkeypatch.setattr(alignment, "entry_points", lambda **kwargs: [type("BadEntry", (), {"load": lambda self: object})()])
    with pytest.raises(ValueError, match="must implement AlignmentEngine"):
        alignment.get_alignment_engine("custom")


@pytest.mark.parametrize("name", ["bowtie2", "bwa-mem3"])
@pytest.mark.parametrize("paired", [False, True])
def test_native_alignment_to_sorted_bam_when_installed(tmp_path, name, paired):
    engine = alignment.get_alignment_engine(name)
    required = {engine.index_command(Path("ref"), 1)[0], name, "samtools"}
    if any(shutil.which(command) is None for command in required):
        pytest.skip("Native alignment tools are not installed")
    sequence = "".join(random.Random(7).choices("ACGT", k=3000))
    ref = tmp_path / "reference.fna"
    ref.write_text(f">chr\n{sequence}\n")
    reads = [tmp_path / "reads 1.fq"]
    if paired:
        reads.append(tmp_path / "reads 2.fq")
    for index, path in enumerate(reads):
        records = []
        for number, start in enumerate((400, 100, 700)):
            read = sequence[start:start + 100] if index == 0 else sequence[start + 150:start + 250].translate(str.maketrans("ACGT", "TGCA"))[::-1]
            records.append(f"@read{number}\n{read}\n+\n{'I' * len(read)}\n")
        path.write_text("".join(records))
    config = WorkflowConfig.legacy(threads=2, sample_count=1)
    runtime = WorkflowRuntime(config, state_dir=tmp_path / "state", logger=WorkflowLogger())
    engine.prepare_index(runtime, ref, sample="native")
    assert engine.index_complete(ref)
    bam = tmp_path / "result.bam"
    runtime.run_shell("alignment", engine.bam_command(ref, reads, 2, bam), sample="native")
    subprocess.run(["samtools", "quickcheck", str(bam)], check=True)
    header = subprocess.check_output(["samtools", "view", "-H", str(bam)], text=True)
    assert "SO:coordinate" in header
    rows = subprocess.check_output(["samtools", "view", str(bam)], text=True).splitlines()
    assert len(rows) == 3 * (2 if paired else 1)
    positions = [int(row.split("\t")[3]) for row in rows]
    assert positions == sorted(positions)
    assert all(not int(row.split("\t")[1]) & 4 for row in rows)
