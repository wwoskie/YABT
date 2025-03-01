import os

from Bio import SeqIO
from Bio import SeqUtils
from numbers import Number
from typing import Iterable


class FastQFilter:
    """
    A utility class for filtering FASTQ files based on quality, GC content, and read length.

    Args:
        path_to_input (str): Path to the input FASTQ file.
        path_to_output (str, optional): Path to the output filtered FASTQ file. If not provided,
            a unique output path will be generated based on the input file name.
        gc_bounds (Iterable[Number] | Number, optional): GC content bounds for filtering.
            Can be a single value or a tuple representing the lower and upper bounds (default: (0, 100)).
        length_bounds (Iterable[Number] | Number, optional): Read length bounds for filtering.
            Can be a single value or a tuple representing the lower and upper bounds (default: (0, 2**32)).
        quality_threshold (Number, optional): Minimum average quality score threshold for filtering
            (default: 0).

    Attributes:
        path_to_input (str): Path to the input FASTQ file.
        path_to_output (str): Path to the output filtered FASTQ file.
        gc_bounds (tuple): Tuple representing the GC content bounds.
        length_bounds (tuple): Tuple representing the read length bounds.
        quality_threshold (Number): Minimum average quality score threshold.

    Methods:
        filter(): Filters the input FASTQ file based on specified criteria.
    """

    def __init__(
        self,
        path_to_input: str,
        path_to_output: str = None,
        gc_bounds: Iterable[Number] | Number = (0, 100),
        length_bounds: Iterable[Number] | Number = (0, 2**32),
        quality_threshold: Number = 0,
    ) -> None:
        self.path_to_input = path_to_input

        if path_to_output is None:
            self.path_to_output = self._uniquify_path(path_to_input)
        else:
            self.path_to_output = path_to_output

        self.gc_bounds = self._make_bounds(gc_bounds)
        self.length_bounds = self._make_bounds(length_bounds)
        self.quality_threshold = quality_threshold

    def filter(self):
        """
        Filters the input FASTQ file based on quality, GC content, and read length.
        Writes the filtered reads to the output FASTQ file.
        """

        reads_that_passed_filtration = []
        for read in SeqIO.parse(self.path_to_input, "fastq"):
            if (
                sum(read.letter_annotations["phred_quality"]) / len(read)
                >= self.quality_threshold
                and self._is_in_bounds(SeqUtils.gc_fraction(read.seq), self.gc_bounds)
                and self._is_in_bounds(len(read), self.length_bounds)
            ):
                reads_that_passed_filtration.append(read)
            SeqIO.write(reads_that_passed_filtration, self.path_to_output, "fastq")

    def _make_bounds(self, bounds) -> tuple:
        """
        Converts bounds to a valid tuple format.

        Args:
            bounds: Bounds (single value or tuple).

        Returns:
            tuple: Bounds in tuple format.
        """

        if not isinstance(bounds, Iterable) and not isinstance(bounds, Number):
            raise TypeError(f"Cannot work with {type(bounds).__name__} type")

        if isinstance(bounds, Number):
            bounds = (0, bounds)
        else:
            bounds = tuple(bounds)

        return bounds

    def _uniquify_path(self, path) -> str:
        """
        Generates a unique output path by appending a counter to the filename.

        Args:
            path (str): Input path.

        Returns:
            str: Unique output path.
        """

        filename, extension = os.path.splitext(path)
        counter = 1

        while os.path.exists(path):
            path = filename + " (" + str(counter) + ")" + extension
            counter += 1

        return path

    def _is_in_bounds(self, value: Number, bounds: tuple) -> bool:
        """
        Checks if a value is within specified bounds.

        Args:
            value (Number): Value to check.
            bounds (tuple): Bounds (lower and upper).

        Returns:
            bool: True if value is within bounds, False otherwise.
        """

        return bounds[0] <= value <= bounds[1]

    def __repr__(self):
        return (
            f"FastQFilter(\n\tpath_to_input={self.path_to_input},"
            + "\n\tpath_to_output={self.path_to_output}"
            + "\n\tgc_bounds={self.gc_bounds},"
            + "\n\tlength_bounds={self.length_bounds},"
            + "\n\tquality_threshold={self.quality_threshold}\n)"
        )
