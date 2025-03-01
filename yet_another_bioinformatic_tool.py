import os

from Bio import SeqIO
from Bio import SeqUtils
from collections import defaultdict
from numbers import Number
from typing import Iterable, Self
from abc import ABC, abstractmethod, abstractproperty


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
            SeqIO.write(reads_that_passed_filtration,
                        self.path_to_output, "fastq")

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


class BiologicalSequence(str, ABC):  # Anton will bully me even more :(
    """
    Abstract base class representing a biological sequence.

    This abstract base class defines common methods and properties for biological sequences.
    Subclasses should implement the `alphabet` property and can override other methods as needed.

    Attributes:
        None

    Properties:
        alphabet (set): The set of characters (symbols) in the biological sequence.

    Methods:
        check_alphabet(alphabet: Iterable) -> bool:
            Checks if the sequence's alphabet is a subset of the given alphabet.
    """

    @property
    @abstractmethod
    def alphabet(self) -> set:
        pass

    def check_alphabet(self, alphabet: Iterable = None) -> bool:
        if alphabet is None:
            alphabet = self.alphabet
        return set(self).issubset(set(alphabet))


class NucleicAcidSequence(BiologicalSequence):
    """
    Represents a nucleic acid sequence (e.g., DNA or RNA).

    This class extends the abstract base class `BiologicalSequence` and provides additional methods
    specific to nucleic acid sequences.

    Attributes:
        None

    Properties:
        alphabet (set): The set of characters (symbols) in the nucleic acid sequence.

    Methods:
        complement() -> str:
            Returns the complement of the nucleic acid sequence using the provided dictionary `_comp_dct`.

        gc_content() -> float:
            Calculates the GC content (percentage of guanine and cytosine bases) in the sequence.

    Raises:
            NotImplementedError: If given instance is of class NucleicAcidSequence.
    """

    @abstractmethod
    def _comp_dct(self):
        """
        Abstract method to define the complementary base pairs dictionary.
        Subclasses must implement this method.

        Returns:
            dict: A dictionary mapping each base to its complement.
        """

        pass

    @property
    def alphabet(self) -> set:  # how can i avoid checking it in every method?
        """
        Returns the set of characters (symbols) in the nucleic acid sequence.

        Returns:
            set: The alphabet of the sequence.
        """

        if not not issubclass(type(self), NucleicAcidSequence):
            raise NotImplementedError(
                "This method is not implemented for NucleicAcidSequence class"
            )
        return set(self._comp_dct)

    def complement(self: Self) -> Self:
        """
        Returns the complement of the nucleic acid sequence using the provided dictionary `_comp_dct`.

        Returns:
            Self: The complement sequence.
        """

        if not issubclass(type(self), NucleicAcidSequence):
            raise NotImplementedError(
                "This method is not implemented for NucleicAcidSequence class"
            )
        return self.__class__("".join([self._comp_dct[letter] for letter in self]))

    @property
    def gc_content(self) -> float:
        """
        Calculates the GC content (percentage of guanine and cytosine bases) in the sequence.

        Returns:
            float: The GC content as a percentage.
        """

        if not not issubclass(type(self), NucleicAcidSequence):
            raise NotImplementedError(
                "This method is not implemented for NucleicAcidSequence class"
            )
        gc_count = 0
        for letter in self:
            if letter.lower() == "g" or letter.lower == "c":
                gc_count += 1
        return gc_count / len(self)


class DNASequence(NucleicAcidSequence):
    """
    Represents a DNA sequence.

    This class extends the `NucleicAcidSequence` and provides methods specific to DNA sequences.

    Attributes:
        None

    Properties:
        _trans_dct (dict): A dictionary mapping DNA bases to their RNA counterparts for transcription.
        _comp_dct (dict): A dictionary mapping DNA bases to their complementary bases.

    Methods:
        transcribe() -> str:
            Transcribes the DNA sequence into RNA by replacing T with U (uracil).
    """

    @property
    def _trans_dct(self) -> dict:
        """
        Returns the dictionary mapping DNA bases to their RNA counterparts for transcription.

        Returns:
            dict: The transcription dictionary.
        """

        return {"T": "U", "u": "t"}

    @property
    def _comp_dct(self) -> dict:
        """
        Returns the dictionary mapping DNA bases to their complementary bases.

        Returns:
            dict: The complement dictionary.
        """
        return {
            "G": "C",
            "C": "G",
            "g": "c",
            "c": "g",
        } | {"A": "T", "T": "A", "a": "t", "t": "a"}

    def transcribe(self) -> str:
        """
        Transcribes the DNA sequence into RNA by replacing T with U (uracil).

        Returns:
            str: The transcribed RNA sequence.
        """

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
    Represents an RNA sequence.

    This class extends the `NucleicAcidSequence` and provides methods specific to RNA sequences.

    Attributes:
        None

    Properties:
        _comp_dct (dict): A dictionary mapping RNA bases to their complementary bases.

    Methods:
        None
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
    Represents an amino acid sequence.

    This class extends the `BiologicalSequence` and provides methods specific to amino acid sequences.

    Attributes:
        None

    Properties:
        alphabet (set): The set of valid amino acid symbols.

    Methods:
        count_aa() -> defaultdict:
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
