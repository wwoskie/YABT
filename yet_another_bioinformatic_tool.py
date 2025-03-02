from pathlib import Path
from typing import Iterable, Iterator

from Bio import SeqIO, SeqUtils
from Bio.SeqRecord import SeqRecord


class FastQFiltrator:
    """
    A class to filter FASTQ reads based on GC content, length, and quality thresholds.

    Parameters
    ----------
    path_to_input : Path | str | None, optional
        Path to the input FASTQ file. If provided, reads will be loaded from this file.
        Default is None.
    reads : Iterable[SeqRecord] | Iterator[SeqRecord] | None, optional
        An iterable or iterator of SeqRecord objects representing the reads.
        If provided, reads will be taken directly from this iterable.
        Default is None.
    gc_bounds : tuple[float, float] | float, optional
        The GC content bounds for filtering. Can be a tuple (lower_bound, upper_bound)
        or a single float representing the upper bound (lower bound defaults to 0).
        Default is (0, 100).
    length_bounds : tuple[int, int] | int, optional
        The length bounds for filtering. Can be a tuple (lower_bound, upper_bound)
        or a single integer representing the upper bound (lower bound defaults to 0).
        Default is (0, 2**32).
    quality_threshold : float, optional
        The minimum average quality score required for a read to pass filtering.
        Default is 0.

    Raises
    ------
    ValueError
        If neither `path_to_input` nor `reads` is provided, or if the provided `reads`
        is not an iterable or iterator.
    TypeError
        If `gc_bounds` or `length_bounds` are not of the expected types.

    Methods
    -------
    filter_fastq()
        Filters the reads based on the specified GC content, length, and quality thresholds.
    write_to_file(path_to_output: Path | str | None = None, rewrite: bool = False)
        Writes the filtered reads to a FASTQ file.
    """

    def __init__(
        self,
        path_to_input: Path | str | None = None,
        reads: Iterable[SeqRecord] | Iterator[SeqRecord] | None = None,
        gc_bounds: tuple[float, float] | float = (0, 100),
        length_bounds: tuple[int, int] | int = (0, 2**32),
        quality_threshold: float = 0,
    ) -> None:
        # Handle read inputs: read from path or take reads
        self.path_to_input = path_to_input
        if isinstance(self.path_to_input, (Path, str)):
            self._read_fastq(self.path_to_input)
        elif isinstance(reads, (Iterable, Iterator)) and self.path_to_input is None:
            self.reads = reads
        else:
            raise ValueError(
                "Reads should be either be: iterable or iterator of reads,"
                + f" but {type(reads).__name__} was given!"
            )

        # Define thresholds
        self.gc_bounds = self._make_bounds(gc_bounds)
        self.length_bounds = self._make_bounds(length_bounds)
        self.quality_threshold = quality_threshold

    def filter_fastq(self) -> None:
        """Filters the reads based on the specified GC content, length, and quality thresholds"""

        reads_that_passed_filtration = []
        for read in self.reads:
            if self._passed_filter(read):
                reads_that_passed_filtration.append(read)
        self.reads = reads_that_passed_filtration

    def write_to_file(
        self, path_to_output: Path | str | None = None, rewrite: bool = False
    ) -> None:
        """
        Writes the filtered reads to a FASTQ file

        Parameters
        ----------
        path_to_output : Path | str | None, optional
            Path to the output FASTQ file. If not provided, the input file path will be used.
            Default is None.
        rewrite : bool, optional
            If True, the output file will be overwritten if it already exists.
            If False, a unique file name will be generated to avoid overwriting.
            Default is False.
        """
        if path_to_output is not None:
            SeqIO.write(
                self.reads, self._uniquify_path(path_to_output, rewrite), "fastq"
            )
        elif path_to_output is None and self.path_to_input is not None:
            path_to_output = self.path_to_input
            SeqIO.write(
                self.reads, self._uniquify_path(path_to_output, rewrite), "fastq"
            )
        else:
            raise ValueError(
                "No output path was provided and could not create it from input path"
            )

    def _passed_filter(self, read: SeqRecord) -> bool:
        return (
            sum(read.letter_annotations["phred_quality"]) / len(read)
            >= self.quality_threshold
            and self._is_in_bounds(SeqUtils.gc_fraction(read.seq), self.gc_bounds)
            and self._is_in_bounds(len(read), self.length_bounds)
        )

    def _read_fastq(self, path_to_input: Path | str) -> None:
        iterator: Iterator[SeqRecord] = SeqIO.parse(path_to_input, "fastq")
        self.reads = iterator

    def _make_bounds(self, bounds: tuple[float, float] | float) -> tuple[float, float]:
        if not (isinstance(bounds, tuple) and len(bounds) == 2) and not isinstance(
            bounds, float
        ):
            raise TypeError(f"Cannot work with {type(bounds).__name__} type")

        elif isinstance(bounds, float):
            bounds = (0, bounds)

        return bounds

    def _uniquify_path(self, path: Path | str, rewrite: bool) -> Path | str:
        if not rewrite:
            path = Path(path) if isinstance(path, str) else path
            parent, filename, extension = path.parent, path.stem, path.suffix
            counter = 1

            while path.exists():
                path = Path(parent, filename + "_(" + str(counter) + ")" + extension)
                counter += 1

        return path

    def _is_in_bounds(self, value: float, bounds: tuple[float, float]) -> bool:
        return bounds[0] <= value <= bounds[1]

    def __str__(self) -> str:
        return (
            "FastQFilter("
            + f"\n\tgc_bounds={self.gc_bounds},"
            + f"\n\tlength_bounds={self.length_bounds},"
            + f"\n\tquality_threshold={self.quality_threshold}\n)"
        )

    def __repr__(self) -> str:
        return (
            f"FastQFilter(\n\tpath_to_input='{self.path_to_input}',"
            + f"\n\tgc_bounds={self.gc_bounds},"
            + f"\n\tlength_bounds={self.length_bounds},"
            + f"\n\tquality_threshold={self.quality_threshold}\n)"
        )
