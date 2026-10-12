import pytest

from metatrawl import matrix_hdf5


@pytest.fixture(params=[0, 1], ids=["zero-based-columns", "one-based-columns"])
def csv_column_base(monkeypatch, request):
    """Simulate CSV default names without changing explicitly named columns."""
    scan_csv = matrix_hdf5.pl.scan_csv

    def scan_with_default_names(*args, **kwargs):
        frame = scan_csv(*args, **kwargs)
        if kwargs.get("new_columns") is None:
            columns = frame.collect_schema().names()
            frame = frame.rename(
                {name: f"column_{request.param + idx}" for idx, name in enumerate(columns)}
            )
        return frame

    monkeypatch.setattr(matrix_hdf5.pl, "scan_csv", scan_with_default_names)


@pytest.mark.parametrize("extra_bed_column", [False, True])
def test_matrix_reference_readers_preserve_coordinates(tmp_path, csv_column_base, extra_bed_column):
    bed_file = tmp_path / "genome.bed"
    rows = ["contigA\t0\t4", "contigA\t4\t10", "contigB\t0\t5", "other\t0\t3"]
    bed_file.write_text("\n".join(row + ("\tinterval" if extra_bed_column else "") for row in rows) + "\n")
    stb_file = tmp_path / "genome.stb"
    stb_file.write_text("contigB\tGENOME\ncontigA\tGENOME\nother\tOTHER\n")
    gene_file = tmp_path / "genes.tsv"
    gene_file.write_text("gene_b\tcontigB\t2\t5\ngene_a\tcontigA\t1\t10\nignored\tother\t1\t3\n")

    scaffolds = matrix_hdf5._collect_scaffold_specs_from_bed_and_stb(
        bed_file=bed_file, stb_file=stb_file, genome="GENOME"
    )
    assert [(s.chrom, s.index_base, s.vector_length, s.min_pos, s.max_pos) for s in scaffolds] == [
        ("contigA", 1, 10, 1, 10),
        ("contigB", 1, 5, 1, 5),
    ]
    genomes, offsets = matrix_hdf5._build_genome_specs(scaffolds)
    assert genomes[0].true_length == 15
    assert genomes[0].matrix_length == 16  # Includes the scaffold separator.
    genes = matrix_hdf5._collect_gene_range_specs(gene_range_table=gene_file, genome_scaffolds=offsets)
    assert [(g.gene, g.chrom, g.axis_start, g.axis_end) for g in genes] == [
        ("gene_a", "contigA", 0, 9),
        ("gene_b", "contigB", 12, 15),
    ]


def test_stb_reader_still_rejects_duplicate_scaffolds(tmp_path, csv_column_base):
    stb_file = tmp_path / "duplicate.stb"
    stb_file.write_text("contigA\tGENOME\ncontigA\tOTHER\n")
    with pytest.raises(ValueError, match="Duplicate scaffold mappings"):
        matrix_hdf5._read_stb_mapping(stb_file)


def test_gene_reader_still_rejects_out_of_bounds_coordinates(tmp_path, csv_column_base):
    bed_file = tmp_path / "genome.bed"
    bed_file.write_text("contigA\t0\t10\n")
    stb_file = tmp_path / "genome.stb"
    stb_file.write_text("contigA\tGENOME\n")
    gene_file = tmp_path / "genes.tsv"
    gene_file.write_text("invalid\tcontigA\t1\t11\n")
    scaffolds = matrix_hdf5._collect_scaffold_specs_from_bed_and_stb(
        bed_file=bed_file, stb_file=stb_file, genome="GENOME"
    )
    _, offsets = matrix_hdf5._build_genome_specs(scaffolds)
    with pytest.raises(ValueError, match="outside the matrix coordinate range"):
        matrix_hdf5._collect_gene_range_specs(gene_range_table=gene_file, genome_scaffolds=offsets)
