import logging

import pytest
from Bio import SeqIO

from yet_another_bioinformatic_tool import FastQFiltrator

logging.disable(logging.CRITICAL)


@pytest.fixture
def example_fastq_path():
    return "example_data/example_fastq.fastq"


@pytest.fixture
def example_reads(example_fastq_path):
    return list(SeqIO.parse(example_fastq_path, "fastq"))


def test_init_with_reads(example_reads):
    filtrator = FastQFiltrator(reads=example_reads, disable_logging=True)
    assert len(filtrator.reads) == 10


def test_init_with_path(example_fastq_path):
    filtrator = FastQFiltrator(path_to_input=example_fastq_path, disable_logging=True)
    assert len(filtrator.reads) == 10


def test_filtering(example_fastq_path):
    filtrator = FastQFiltrator(
        path_to_input=example_fastq_path,
        gc_bounds=(0, 60),
        length_bounds=(80, 300),
        disable_logging=True,
    )

    filtrator.filter_fastq()
    assert len(filtrator.reads) == 8


def test_init_with_both(example_reads, example_fastq_path):
    with pytest.raises(ValueError):
        FastQFiltrator(
            path_to_input=example_fastq_path, reads=example_reads, disable_logging=True
        )


def test_invalid_gc_bounds_type(example_reads):
    with pytest.raises(TypeError):
        FastQFiltrator(reads=example_reads, gc_bounds="some str", disable_logging=True)


def test_invalid_length_bounds_type(example_reads):
    with pytest.raises(TypeError):
        FastQFiltrator(
            reads=example_reads, length_bounds="some str", disable_logging=True
        )


def test_write_to_file(tmp_path, example_reads):
    output_path = tmp_path / "output.fastq"
    filtrator = FastQFiltrator(reads=example_reads, disable_logging=True)
    filtrator.filter_fastq()
    filtrator.write_to_file(output_path)
    assert output_path.exists()
    assert len(list(SeqIO.parse(output_path, "fastq"))) == 10
