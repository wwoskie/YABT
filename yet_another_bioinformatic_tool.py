from abc import ABC, abstractmethod
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Iterable, Iterator, Self

from Bio import SeqIO, SeqUtils
from Bio.SeqRecord import SeqRecord

from logger import setup_class_logger


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
    gc_bounds : Iterable[float, float] | float, optional
        The GC content bounds for filtering. Can be a tuple (lower_bound, upper_bound)
        or a single float representing the upper bound (lower bound defaults to 0).
        Default is (0, 1).
    length_bounds : Iterable[int, int] | int, optional
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
        reads: list[SeqRecord]
        | Iterable[SeqRecord]
        | Iterator[SeqRecord]
        | None = None,
        gc_bounds: tuple[float, float] | float = (0, 1),
        length_bounds: tuple[int, int] | int = (0, 2**32),
        quality_threshold: float = 0,
        logs_dir: Path | str = "",
        log_file: Path | str = "auto",
        disable_logging: bool = False,
    ) -> None:
        # Handle read inputs: read from path or take reads
        # Setup logging
        if log_file == "auto":
            log_file = f"{datetime.now().strftime('%Y%m%d-%H%M%S-%f')[:-3]}-{self.__class__.__name__}.log"

        if logs_dir is not None and log_file is not None:
            log_file = Path(
                logs_dir,
                log_file,
            )
        else:
            log_file = None

        if disable_logging:
            log_file = None

        self.logger = setup_class_logger(
            self.__class__.__name__,
            verbosity_console=2,
            verbosiy_file=2,
            log_file=log_file,
        )

        self.path_to_input = path_to_input
        if isinstance(self.path_to_input, (Path, str)) and (
            isinstance(reads, list) or isinstance(reads, (Iterable, Iterator))
        ):
            raise ValueError(
                "Reads should either be path or reads iterator/list, not both"
            )
        elif isinstance(self.path_to_input, (Path, str)):
            self.reads = self._read_fastq(self.path_to_input)
        elif isinstance(reads, (Iterable, Iterator)) and self.path_to_input is None:
            self.reads = list(reads)
        elif isinstance(reads, list) and self.path_to_input is None:
            self.reads = reads
        else:
            raise ValueError(
                "Reads should be either be: list of reads,"
                + f" but {type(reads).__name__} was given!"
            )

        # Define thresholds
        try:
            self.gc_bounds = self._make_bounds(gc_bounds)
        except TypeError as e:
            self.logger.error(e, stack_info=True, exc_info=True)
            raise e

        try:
            self.length_bounds = self._make_bounds(length_bounds)
        except TypeError as e:
            self.logger.error(e, stack_info=True, exc_info=True)
            raise e

        self.quality_threshold = quality_threshold

        self.logger.info(f"Starting fastq filtering: {self}")

    def filter_fastq(self) -> None:
        """Filters the reads based on the specified GC content, length, and quality thresholds"""
        self.logger.info("Filtering fastq...")
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
        if isinstance(path_to_output, (str, Path)):
            path_to_output = path_to_output
        elif path_to_output is None and self.path_to_input is not None:
            path_to_output = self.path_to_input
        else:
            raise ValueError(
                "No output path was provided and could not create it from input path"
            )

        self.logger.info(f"Writing fastq to file: {path_to_output}...")

        SeqIO.write(self.reads, self._uniquify_path(path_to_output, rewrite), "fastq")

    def _passed_filter(self, read: SeqRecord) -> bool:
        return (
            sum(read.letter_annotations["phred_quality"]) / len(read)
            >= self.quality_threshold
            and self._is_in_bounds(SeqUtils.gc_fraction(read.seq), self.gc_bounds)
            and self._is_in_bounds(len(read), self.length_bounds)
        )

    def _read_fastq(self, path_to_input: Path | str) -> list:
        iterator: Iterator[SeqRecord] = list(SeqIO.parse(path_to_input, "fastq"))
        return iterator

    def _make_bounds(self, bounds: tuple[float, float] | float) -> tuple[float, float]:
        if not (isinstance(bounds, Iterable) and len(bounds) == 2) and not isinstance(
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
            f"FastQFiltrator(path_to_input='{self.path_to_input}',"
            + f" gc_bounds={self.gc_bounds},"
            + f" length_bounds={self.length_bounds},"
            + f" quality_threshold={self.quality_threshold})"
        )

    def __repr__(self) -> str:
        return (
            f"FastQFiltrator(path_to_input='{self.path_to_input}',"
            + f" gc_bounds={self.gc_bounds},"
            + f" length_bounds={self.length_bounds},"
            + f" quality_threshold={self.quality_threshold})"
        )


class BiologicalSequence(str, ABC):
    """
    Abstract base class representing a biological sequence.

    Attributes
    ----------
    alphabet : set
        The set of valid characters in the sequence.

    Methods
    -------
    check_alphabet(alphabet: Iterable = None) -> bool
        Checks if the sequence characters are within the given alphabet.
    """

    @property
    @abstractmethod
    def alphabet(self) -> set[str]:
        pass

    def check_alphabet(self, alphabet: Iterable[str] | None = None) -> bool:
        if alphabet is None:
            alphabet = self.alphabet
        return set(self).issubset(set(alphabet))


class NucleicAcidSequence(BiologicalSequence):
    """
    Abstract base class representing a nucleic acid sequence.

    Methods
    -------
    complement() -> Self
        Returns the complement of the sequence.
    gc_content() -> float
        Calculates the GC content of the sequence.
    """

    @abstractmethod
    def _comp_dct(self) -> None:
        self._check_implementation()
        pass

    @property
    def alphabet(self) -> set[str]:
        self._check_implementation()
        return set(self._comp_dct)

    def complement(self: Self) -> Self:
        self._check_implementation()
        return self.__class__("".join([self._comp_dct[letter] for letter in self]))

    @property
    def gc_content(self) -> float:
        self._check_implementation()
        gc_count = 0
        for letter in self:
            if letter.lower() == "g" or letter.lower() == "c":
                gc_count += 1
        return gc_count / len(self)

    def _check_implementation(self) -> None:
        if type(self).__name__ == "NucleicAcidSequence":
            raise NotImplementedError(
                "This method is not implemented for 'NucleicAcidSequence'"
            )


class DNASequence(NucleicAcidSequence):
    """
    Class representing a DNA sequence.

    Methods
    -------
    transcribe() -> RNASequence
        Transcribes the DNA sequence into an RNA sequence.
    """

    @property
    def _trans_dct(self) -> dict:
        return {"T": "U", "u": "t"}

    @property
    def _comp_dct(self) -> dict:
        return {
            "G": "C",
            "C": "G",
            "g": "c",
            "c": "g",
        } | {"A": "T", "T": "A", "a": "t", "t": "a"}

    def transcribe(self) -> str:
        return RNASequence(
            "".join(
                [
                    self._trans_dct[letter] if letter in self._trans_dct else letter
                    for letter in self
                ]
            )
        )

    def __repr__(self):
        return f"DNASequence('{self}')"


class RNASequence(NucleicAcidSequence):
    """
    Class representing an RNA sequence.
    """

    @property
    def _comp_dct(self) -> dict:
        return {
            "G": "C",
            "C": "G",
            "g": "c",
            "c": "g",
        } | {"A": "U", "U": "A", "a": "u", "u": "a"}

    def __repr__(self):
        return f"RNASequence('{self}')"


class AminoAcidSequence(BiologicalSequence):
    """
    Class representing an amino acid sequence.

    Methods
    -------
    count_aa() -> defaultdict
        Counts the occurrences of each amino acid in the sequence.
    """

    __alphabet_cap = "FLSYCWPHQRIMTNKVADEGUO"
    _alphabet = __alphabet_cap + __alphabet_cap.lower()

    @property
    def alphabet(self) -> set:
        return set(self._alphabet)

    def count_aa(self) -> defaultdict:
        aa_counts = defaultdict(int)

        for aa in self:
            aa_counts[aa] += 1

        return aa_counts

    def __repr__(self):
        return f"AminoAcidSequence('{self}')"
