import os
import functools
import sys
import re
import requests
import warnings

from abc import ABC, abstractmethod, abstractproperty
from Bio import SeqIO, SeqUtils
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timedelta
from io import StringIO
from numbers import Number
from typing import Callable, Dict, Iterable, List, Self

TG_API_TOKEN = os.getenv("TG_API_TOKEN")
CHAT_ID = os.getenv("CHAT_ID")


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


def run_genscan(
    sequence: str = None,
    sequence_file: str = None,
    organism: str = "Vertebrate",
    exon_cutoff: float = 1.00,
    sequence_name: str = "",
) -> GenscanOutput:
    """
    API to access GENSCAN online tool

    Args:
        sequence (str): The input DNA sequence to be analyzed by GENSCAN.
        sequence_file (str): The path to a file containing the input DNA sequence in fasta format.
        organism (str): The type of organism for which gene prediction is performed. Available organisms are: 'Arabidopsis', 'Maize', 'Vertebrate' (default is 'Vertebrate').
        exon_cutoff (float): The cutoff value for exons predicted by GENSCAN (default is 1.00).
        sequence_name (str): The name of the input sequence.

    Returns:
        GenscanOutput: An instance of the GenscanOutput class containing the output data from the GENSCAN prediction.
    """

    url = "http://argonaute.mit.edu/cgi-bin/genscanw_py.cgi"

    if organism not in {"Arabidopsis", "Maize", "Vertebrate"}:
        raise ValueError(
            "Organism not available. Choose one from the list: 'Arabidopsis', 'Maize', 'Vertebrate'"
        )

    if sequence_file is not None:
        with open(f"{sequence_file}") as f:
            file = f.read()
    else:
        file = None

    if sequence is not None and sequence_file is not None:
        warnings.warn(
            "You've passed both sequence file and sequence string, sequence string will be used",
            UserWarning,
        )

    data = {
        "-o": organism,
        "-e": exon_cutoff,
        "-n": sequence_name,
        "-p": "Predicted peptides only",
        "-u": file,
        "-s": sequence,
    }

    p = requests.post(url=url, data=data)

    output = p.text.split("<pre>")[1].split("</pre>")[0]

    output_list = list(
        filter(
            None,
            re.split(
                "fasta : | bp|Predicted genes/exons:"
                + "|Suboptimal exons with probability|"
                + f"Exnum|Predicted peptide sequence{re.escape('(s):')}",
                output,
            ),
        )
    )

    status = p.status_code

    exon_dict = process_table(output_list[3].strip().split("\n\n\n\n")[1:])

    suboptimal_exon_dict = process_table(output_list[5].strip().split("\n\n\n\n")[1:])

    exon_dict |= suboptimal_exon_dict

    intron_dict = compute_introns(exon_dict)

    cds_list = process_proteins(output_list[6].strip())

    output_dict = {
        "status": status,
        "sequence_name": sequence_name,
        "exon_dict": exon_dict,
        "intron_dict": intron_dict,
        "cds_dict": cds_list,
    }

    return GenscanOutput(**output_dict)


def telegram_logger(chat_id: int) -> Callable:
    """
    Decorator that logs the execution of another function and sends the log (if available) to a specified Telegram chat. It takes a chat_id as input and returns a Callable.

    Args:
        chat_id (int): The chat id of the Telegram user where the log will be sent.

    Returns:
        Callable: A decorator function that can be used to log the execution of other functions.
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):

            output = None

            try:
                start_time = datetime.now()

                with Capturing() as captured:
                    output = func(*args, **kwargs)

                end_time = datetime.now()

                delta = end_time - start_time
                # delta = delta + timedelta(days=10000)

                text = (
                    f"😎👍💅🏼\nFunction `{func.__name__}` "
                    + f"successfully executed in `{delta}`"
                )

            except Exception as e:
                text = (
                    f"😳😢🤔\nFunction `{func.__name__}` "
                    + f"failed with an exception:\n\n `{type(e).__name__}: {e}`"
                )

            finally:
                token = TG_API_TOKEN

                base_url = "https://api.telegram.org/bot" + token

                if len(captured.getvalue()) > 0:

                    p = requests.post(
                        url=base_url + "/sendDocument",
                        data={"chat_id": chat_id},
                        files={
                            "document": (f"{func.__name__}.log", captured.getvalue())
                        },
                    )

                r = requests.get(
                    url=base_url + "/sendMessage",
                    params={"chat_id": chat_id, "text": text, "parse_mode": "markdown"},
                )

            return output

        return wrapper

    return decorator
