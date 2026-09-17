"""Database compatibility checks using small real DuckDB/Parquet fixtures."""
from dataclasses import replace
import json
from pathlib import Path
import shutil
import uuid

from click.testing import CliRunner
import duckdb
import polars as pl
import pytest

from metatrawl import db, provenance, workflows, cache
from metatrawl.cli import cli
from metatrawl.config import ProfileConfig, WorkflowConfig


def bundle(root, sample, *, config=None, manifest=True, references=None, details=None):
    profile = root / f"{sample}.profile.parquet"
    stats = root / f"{sample}.genome_stats.parquet"
    pl.DataFrame(dict(chrom=["c"], pos=[1], genome=["g"], A=[5], C=[0], G=[0], T=[0])).write_parquet(profile)
    pl.DataFrame(dict(genome=["g"], coverage=[5.0], breadth=[1.0], ber=[1.0])).write_parquet(stats)
    if manifest:
        provenance.write_manifest(profile_file=profile, sample_id=sample,
                                  contract=provenance.profiling_contract(config or ProfileConfig()),
                                  references=references or {}, details=details or {})
    return db.ProfileBundle(sample, profile, stats)


@pytest.fixture
def conn(tmp_path):
    with db.connect(tmp_path / "project.duckdb") as connection:
        db.add_runs(connection, ["S1", "S2", "S3"])
        yield connection


def test_uuid_survives_reopen_move_copy_and_new_db_differs(tmp_path):
    path = tmp_path / 'a.duckdb'
    with db.connect(path) as c:
        first = provenance.describe(c)
        uuid.UUID(first['database_uuid'])
    copy = tmp_path / 'copy.duckdb'
    shutil.copy2(path, copy)
    renamed = path.rename(tmp_path / 'renamed.duckdb')
    for p in [renamed, copy]:
        with db.connect(p) as c:
            assert provenance.describe(c)['database_uuid'] == first['database_uuid']
    with db.connect(path) as c:
        assert provenance.describe(c)['database_uuid'] != first['database_uuid']


def test_read_only_legacy_does_not_migrate_and_adoption_preserves_payload_types(tmp_path):
    path = tmp_path / 'legacy.duckdb'
    with duckdb.connect(str(path)) as c:
        c.execute(db.SCHEMA_SQL)
        c.execute("INSERT INTO samples VALUES ('S1','S1','complete',1,1)")
        c.execute("ALTER TABLE profile_positions ALTER COLUMN A TYPE DOUBLE")
        c.execute("INSERT INTO profile_positions VALUES ('S1','c',1,'g',5,0,0,0,NULL)")
    with db.connect_read_only(path) as c:
        assert provenance.describe(c)['schema_version'] == 'legacy-unversioned'
        assert not provenance.table_exists(c, 'database_metadata')
    with db.connect(path) as c:
        assert provenance.describe(c)['origin'] == 'legacy-adopted'
        assert c.execute("SELECT provenance_status FROM sample_provenance").fetchone()[0] == 'legacy-unknown'
        assert c.execute("SELECT typeof(A),A FROM profile_positions").fetchone() == ('DOUBLE', 5.0)
        identity = provenance.describe(c)['database_uuid']
    with db.connect(path) as c:
        assert provenance.describe(c)['database_uuid'] == identity
        assert c.execute("SELECT count(*) FROM sample_provenance").fetchone()[0] == 1


@pytest.mark.parametrize('read_only', [False, True])
def test_future_schema_refused_before_ddl(tmp_path, read_only):
    path = tmp_path / 'newer.duckdb'
    with db.connect(path) as c:
        c.execute('UPDATE database_metadata SET schema_version=999')
        c.execute('DROP TABLE gene_stats')
    with pytest.raises(ValueError, match='Unsupported.*schema'):
        (db.connect_read_only if read_only else db.connect)(path)
    with duckdb.connect(str(path), read_only=True) as c:
        assert not provenance.table_exists(c, 'gene_stats')


def test_identical_contract_import_and_no_global_profile_scan(conn, tmp_path):
    statements=[]
    class Recorder:
        def execute(self, sql, *args):
            statements.append(sql)
            return conn.execute(sql, *args)
    for name in ['S1','S2']:
        db.import_profile_bundle(Recorder(), bundle(tmp_path, name))
    assert conn.execute('SELECT count(*) FROM profiling_contracts').fetchone()[0] == 1
    assert conn.execute('SELECT count(*) FROM sample_provenance').fetchone()[0] == 2
    assert not any('FROM PROFILE_POSITIONS' in sql.upper() for sql in statements)


def test_mismatch_rejected_before_payload_insert_and_override_is_audited(conn, tmp_path):
    db.import_profile_bundle(conn, bundle(tmp_path,'S1'))
    original = provenance.describe(conn)['accepted_contract_id']
    b = bundle(tmp_path,'S2',config=replace(ProfileConfig(),min_freq=0.1))
    with pytest.raises(ValueError,match='profile.min_freq'):
        db.import_profile_bundle(conn,b)
    assert conn.execute("SELECT count(*) FROM samples WHERE sample_id='S2'").fetchone()[0] == 0
    db.import_profile_bundle(conn,b,allow_incompatible_profiles=True)
    assert provenance.describe(conn)['accepted_contract_id'] == original
    assert conn.execute("SELECT count(*) FROM provenance_events WHERE kind='compatibility-override'").fetchone()[0] == 1
    assert conn.execute('SELECT count(*) FROM profiling_contracts').fetchone()[0] == 2


def test_unknown_import_needs_override_and_is_never_relabelled(conn,tmp_path):
    b=bundle(tmp_path,'S1',manifest=False)
    with pytest.raises(ValueError,match='unknown historical'):
        db.import_profile_bundle(conn,b)
    db.import_profile_bundle(conn,b,allow_incompatible_profiles=True)
    assert conn.execute('SELECT contract_id,provenance_status FROM sample_provenance').fetchone() == (None,'legacy-unknown')
    with pytest.raises(ValueError,match='legacy-unknown'):
        db.import_profile_bundle(conn,bundle(tmp_path,'S2'))
    db.import_profile_bundle(conn,bundle(tmp_path,'S2'),allow_incompatible_profiles=True)
    assert conn.execute("SELECT contract_id FROM sample_provenance WHERE sample_id='S1'").fetchone()[0] is None


def test_legacy_declaration_does_not_claim_verification(conn,tmp_path):
    db.import_profile_bundle(conn,bundle(tmp_path,'S1',manifest=False),allow_incompatible_profiles=True)
    c=provenance.profiling_contract(ProfileConfig())
    assert provenance.declare_legacy(conn,contract=c,samples=['S1'],reason='Archived configuration') == 1
    assert conn.execute('SELECT provenance_status FROM sample_provenance').fetchone()[0] == 'user-declared'
    db.import_profile_bundle(conn,bundle(tmp_path,'S2'))
    with pytest.raises(ValueError,match='legacy-unknown'):
        provenance.declare_legacy(conn,contract=c,samples=['S2'],reason='invalid')


def test_payload_failure_rolls_back_all_provenance(conn,tmp_path):
    a=bundle(tmp_path,'S1'); b=bundle(tmp_path,'S2')
    b.genome_stats_file.write_bytes(b'broken parquet')
    with pytest.raises(duckdb.Error):
        db.import_profile_bundles(conn,[a,b])
    for table in ['profile_positions','samples','sample_provenance','profiling_contracts','profile_references']:
        assert conn.execute(f'SELECT count(*) FROM {table}').fetchone()[0] == 0
    assert provenance.describe(conn)['accepted_contract_id'] is None


@pytest.mark.parametrize('key', ['sample_id','contract_id','manifest_version','status'])
def test_malformed_manifest_cannot_be_forced(conn,tmp_path,key):
    b=bundle(tmp_path,'S1')
    path=provenance.manifest_path(b.profile_file)
    value=json.loads(path.read_text()); value[key]='invalid'; path.write_text(json.dumps(value))
    with pytest.raises(ValueError):
        db.import_profile_bundle(conn,b,allow_incompatible_profiles=True)


def test_explicit_missing_manifest_is_not_silently_unknown(conn,tmp_path):
    b=replace(bundle(tmp_path,'S1'),provenance_file=tmp_path/'absent.json')
    with pytest.raises(FileNotFoundError):
        db.import_profile_bundle(conn,b,allow_incompatible_profiles=True)


def test_reference_difference_is_checked(conn,tmp_path):
    reference={'sequence_hash':'abc','annotation_hash':'gene1','scaffolds':[['c',5,'abc']]}
    db.import_profile_bundle(conn,bundle(tmp_path,'S1',references={'g':reference}))
    changed={**reference,'annotation_hash':'gene2'}
    with pytest.raises(ValueError,match='annotation differ'):
        db.import_profile_bundle(conn,bundle(tmp_path,'S2',references={'g':changed}))
    db.import_profile_bundle(conn,bundle(tmp_path,'S2',references={'g':changed}),allow_incompatible_profiles=True)
    assert conn.execute('SELECT count(*) FROM profile_references').fetchone()[0] == 2


def test_null_model_content_difference_is_checked(conn,tmp_path):
    db.import_profile_bundle(conn,bundle(tmp_path,'S1',details={'null_model_sha256':'one'}))
    with pytest.raises(ValueError,match='null-model content'):
        db.import_profile_bundle(conn,bundle(tmp_path,'S2',details={'null_model_sha256':'two'}))


def test_preflight_rejects_before_work_and_keeps_baseline(conn):
    c=provenance.profiling_contract(ProfileConfig())
    run=provenance.start_run(conn,contract=c,details={'threads':1})
    provenance.finish_run(conn,run,'complete')
    with pytest.raises(ValueError,match='min_baseq'):
        provenance.start_run(conn,contract=provenance.profiling_contract(replace(ProfileConfig(),min_baseq=30)),details={})
    assert conn.execute('SELECT count(*) FROM workflow_runs').fetchone()[0] == 1
    assert conn.execute('SELECT status FROM workflow_runs').fetchone()[0] == 'complete'


def test_operational_changes_not_scientific_and_environment_redacted():
    config=WorkflowConfig.legacy(threads=2,sample_count=2)
    config2=replace(config,sample_workers=50)
    assert provenance.profiling_contract(config.profile) == provenance.profiling_contract(config2.profile)
    name=next(iter(config.stages))
    config.stages[name]=replace(config.stages[name],environment={'TOKEN':'secret'})
    assert 'secret' not in provenance.canonical(provenance.operational_config(config))


def test_references_ignore_concatenation_order_but_detect_sequence_or_gene_changes(tmp_path):
    ref=cache.PreparedReference(tmp_path/'ref.fa',tmp_path/'genes.fa',tmp_path/'ref.stb')
    ref.reference_fasta.write_text('>b\nACGT\n>a\nTT\n')
    ref.gene_fasta.write_text('>b_1 # 1 # 2 # 1\nAC\n>a_1 # 1 # 2 # 1\nTT\n')
    ref.stb_file.write_text('a\tg1\nb\tg2\n')
    first=provenance.reference_catalog(ref)
    ref.reference_fasta.write_text('>a description\ntt\n>b\nAC\nGT\n')
    assert provenance.reference_catalog(ref) == first
    ref.gene_fasta.write_text('>b_1 # 1 # 2 # 1\nAG\n>a_1 # 1 # 2 # 1\nTT\n')
    second=provenance.reference_catalog(ref)
    assert second['g1']==first['g1']
    assert second['g2']['annotation_hash']!=first['g2']['annotation_hash']


def test_cleanup_deletes_manifest_only_after_success(conn,tmp_path):
    b=bundle(tmp_path,'S1')
    db.import_profile_bundle(conn,b)
    assert workflows.cleanup_profile_bundle(b) == 3
    assert not provenance.manifest_path(b.profile_file).exists()
    assert json.loads(conn.execute('SELECT manifest_json FROM sample_provenance').fetchone()[0])['sample_id']=='S1'


def test_database_info_cli_is_read_only(tmp_path):
    path=tmp_path/'db.duckdb'
    with duckdb.connect(str(path)) as c: c.execute(db.SCHEMA_SQL)
    result=CliRunner().invoke(cli,['database','info','--db',str(path)])
    assert result.exit_code==0, result.output
    assert 'legacy-unversioned' in result.output
    with duckdb.connect(str(path)) as c: assert not provenance.table_exists(c,'database_metadata')


def test_sync_preflight_refuses_before_workers(tmp_path,monkeypatch):
    path=tmp_path/'db.duckdb'
    with db.connect(path) as c:
        db.add_runs(c,['S1','S2'])
        db.import_profile_bundle(c,bundle(tmp_path,'S1'))
    def forbidden(**kwargs): raise AssertionError('Worker must not start')
    monkeypatch.setattr(workflows,'profile_sra_runs',forbidden)
    config=WorkflowConfig.legacy(threads=1,sample_count=1)
    config=replace(config,profile=replace(config.profile,min_freq=.5))
    with pytest.raises(ValueError,match='min_freq'):
        workflows.sync_remaining_profiles(db_file=path,cache_dir=tmp_path/'cache',scratch_dir=tmp_path/'scratch',output_dir=tmp_path/'outputs',workflow_config=config,check_dependencies=False)


def test_allele_reference_cache_is_cleared_on_batch_rollback(tmp_path):
    from test_allele_mask import _project_files, _write_bundle
    from metatrawl import allele_mask
    cache_dir, _, _, _ = _project_files(tmp_path)
    a = _write_bundle(tmp_path, 'S1')
    b = _write_bundle(tmp_path, 'S2')
    b.genome_stats_file.write_bytes(b'broken parquet')
    reference_cache = allele_mask.AlleleMaskWriteCache()
    with db.connect(tmp_path/'mask.duckdb') as c:
        db.configure_profile_storage(c, mode='allele-mask', min_cov=5)
        db.add_runs(c, ['S1','S2'])
        with pytest.raises(duckdb.Error):
            db.import_profile_bundles(c, [a,b], cache_dir=cache_dir, reference_cache=reference_cache)
        assert not reference_cache._ensured_genomes
        assert not reference_cache._segments
        db.import_profile_bundle(c, a, cache_dir=cache_dir, reference_cache=reference_cache)
        assert c.execute('SELECT count(*) FROM allele_mask_reference_segments').fetchone()[0] > 0
        assert c.execute('SELECT count(*) FROM allele_mask_profile_blocks b LEFT JOIN allele_mask_reference_segments r USING(segment_id) WHERE r.segment_id IS NULL').fetchone()[0] == 0


def test_batch_contract_mismatch_rolls_back_earlier_samples(conn,tmp_path):
    a=bundle(tmp_path,'S1'); b=bundle(tmp_path,'S2',config=replace(ProfileConfig(),min_freq=.1))
    with pytest.raises(ValueError,match='min_freq'):
        db.import_profile_bundles(conn,[a,b])
    assert conn.execute('SELECT count(*) FROM profile_positions').fetchone()[0]==0
    assert conn.execute('SELECT count(*) FROM sample_provenance').fetchone()[0]==0


def test_checkpoint_mismatch_rejects_even_with_import_override(tmp_path):
    scratch=tmp_path/'scratch'; scratch.mkdir()
    previous=provenance.profiling_contract(ProfileConfig())
    (scratch/'profiling_settings.json').write_text(json.dumps(previous))
    ref=cache.PreparedReference(tmp_path/'reference.fa',tmp_path/'genes.fa',tmp_path/'reference.stb')
    with pytest.raises(ValueError,match='checkpoint settings differ'):
        workflows._run_alignment_and_profile(run_id='S1',sample_scratch=scratch,reference=ref,
            output_dir=tmp_path/'out',profile_config=replace(ProfileConfig(),min_freq=.8),
            logger=workflows.WorkflowLogger(),allow_incompatible_profiles=True)


def test_old_bam_checkpoint_requires_acknowledgement(tmp_path):
    scratch=tmp_path/'scratch'; scratch.mkdir()
    (scratch/'S1.bam').write_bytes(b'old')
    ref=cache.PreparedReference(tmp_path/'reference.fa',tmp_path/'genes.fa',tmp_path/'reference.stb')
    with pytest.raises(ValueError,match='unknown historical settings'):
        workflows._run_alignment_and_profile(run_id='S1',sample_scratch=scratch,reference=ref,
            output_dir=tmp_path/'out',logger=workflows.WorkflowLogger())
    assert not (scratch/'profiling_settings.json').exists()


def test_database_cli_template_declaration_and_history(tmp_path):
    path=tmp_path/'project.duckdb'
    with db.connect(path) as c:
        db.add_runs(c,['S1'])
        db.import_profile_bundle(c,bundle(tmp_path,'S1',manifest=False),allow_incompatible_profiles=True)
    runner=CliRunner()
    template=runner.invoke(cli,['database','contract-template'])
    assert template.exit_code==0,template.output
    contract_file=tmp_path/'contract.json'; contract_file.write_text(template.output)
    result=runner.invoke(cli,['database','declare-contract','--db',str(path),
                             '--contract-file',str(contract_file),'--reason','Archived config'])
    assert result.exit_code==0,result.output
    assert 'user-declared' in result.output
    history=runner.invoke(cli,['database','history','--db',str(path)])
    assert history.exit_code==0
    assert 'legacy-contract-declared' in history.output


def test_user_declared_samples_do_not_become_recorded_on_reopen(conn,tmp_path):
    db.import_profile_bundle(conn,bundle(tmp_path,'S1',manifest=False),allow_incompatible_profiles=True)
    provenance.declare_legacy(conn,contract=provenance.profiling_contract(ProfileConfig()),samples=None,reason='Config')
    db.init_schema(conn)
    assert provenance.describe(conn)['sample_provenance']=={'user-declared':1}


def test_schema_migration_never_selects_from_profile_positions(tmp_path):
    with duckdb.connect(str(tmp_path/'legacy.duckdb')) as c:
        c.execute(db.SCHEMA_SQL)
        statements=[]
        class Recorder:
            def execute(self,sql,*args):
                statements.append(sql)
                return c.execute(sql,*args)
        db.init_schema(Recorder())
        assert not any('FROM PROFILE_POSITIONS' in sql.upper() for sql in statements)
        assert not any('ALTER COLUMN' in sql.upper() and 'PROFILE_POSITIONS' in sql.upper() for sql in statements)


def test_storage_conversion_preserves_provenance_with_derived_identity(tmp_path):
    from test_allele_mask import _project_files, _make_database
    from metatrawl import migration
    cache_dir, _, _, _ = _project_files(tmp_path)
    source = _make_database(tmp_path, name='source', mode='full', cache_dir=cache_dir)
    target = tmp_path/'converted.duckdb'
    with db.connect_read_only(source) as c:
        source_metadata=provenance.describe(c)
        samples=c.execute('SELECT sample_id, contract_id, provenance_status FROM sample_provenance ORDER BY sample_id').fetchall()
    migration.migrate_full_database(source_db=source,output_db=target,cache_dir=cache_dir,min_cov=5)
    with db.connect_read_only(target) as c:
        metadata=provenance.describe(c)
        assert metadata['database_uuid']!=source_metadata['database_uuid']
        assert metadata['accepted_contract_id']==source_metadata['accepted_contract_id']
        assert c.execute('SELECT sample_id, contract_id, provenance_status FROM sample_provenance ORDER BY sample_id').fetchall()==samples
        assert c.execute("SELECT count(*) FROM provenance_events WHERE kind='storage-conversion'").fetchone()[0]==1


def test_invalid_contract_parameters_cannot_be_forced(conn,tmp_path):
    b=bundle(tmp_path,'S1')
    path=provenance.manifest_path(b.profile_file)
    value=json.loads(path.read_text())
    value['contract']['profile']['min_freq']=2.0
    value['contract_id']=provenance.fingerprint(value['contract'])
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError,match='min_freq'):
        db.import_profile_bundle(conn,b,allow_incompatible_profiles=True)
