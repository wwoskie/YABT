import os

from Bio import SeqIO
from Bio import SeqUtils
from numbers import Number
from collections.abc import Iterable


import modules.dna_rna_tools as dna_rna_tools
import modules.protein_tools as protein_tools
import modules.fastq_tools as fastq_tools
from modules.dna_rna_tools import create_input_dict
from modules.protein_tools import read_seq_from_fasta


command_dict_nucl = {
    "check_seq_type": dna_rna_tools.check_seq_type,
    "reverse": dna_rna_tools.reverse,
    "complement": dna_rna_tools.complement,
    "transcribe": dna_rna_tools.transcribe,
}


command_dct_prot = {
    "find_sites": protein_tools.find_sites,
    "get_protein_rnas_number": protein_tools.get_protein_rnas_number,
    "is_protein_valid": protein_tools.is_protein_valid,
    "get_length_of_protein": protein_tools.get_length_of_protein,
    "count_aa": protein_tools.count_aa,
    "get_fracture_of_aa": protein_tools.get_fracture_of_aa,
}


def run_dna_rna_tools(seqs: dict, command: str) -> dict:
    """
    Runs dna_rna_tools on given dict of seqs with given command

    Arguments:
    - seqs (dict): Input dict of format {seq_name: 'seq'}

    Return:
    - output_dict (dict | str | bool):
        dict of results of operations {seq_name: 'result'} or one operation result
        if one sequence given
    """

    output_dict = {}

    for seq_name, seq in seqs.items():
        if command == "check_seq_type":  # user may want to check given seqs
            output_dict |= {seq_name: dna_rna_tools.check_seq_type(seq)}
        else:  # if other command
            nucl_type = dna_rna_tools.check_seq_type(seq)
            if nucl_type is None:
                raise ValueError("Can only work with DNA or RNA sequence")

            if command == "complement":
                output_dict |= {seq_name: dna_rna_tools.complement(seq, nucl_type)}
            else:
                output_dict |= {seq_name: command_dict_nucl[command](seq)}

    if len(output_dict) == 1:
        return output_dict[list(output_dict.keys())[0]]
    return output_dict


def run_fastq_tools(
    seqs: dict,  # how can i make native type hint here?
    gc_bounds: tuple | list | int | float = (0, 100),
    length_bounds: tuple | list | int | float = (0, 2**32),
    quality_threshold: int | float = 0,
) -> dict:
    """
    Runs fastq filtration by GC-content, length and quality procedure on input
    dict of format {'seq_name': ('nucl_seq', 'quality_for_seq')}

    Arguments:
    - seqs (dict): input dict of format {'seq_name': ('nucl_seq', 'quality_for_seq')}
    - gc_bounds (tuple | list | int | float):
        bounds for GC filtration, can process int, float, or tuple and list of
        length 2. Default is (0, 100) (not filtered by GC)
    - length_bounds (tuple | list| int | float):
        bounds for length filtration, full analog of gc_bounds.
        Default is (0, 2**32)
    - quality_threshold (int | float):
        quality threshold to check against. Default is 0

    Return:
    - passed_filtration_seqs (dict):
        dict of format {'seq_name': ('nucl_seq', 'quality_for_seq')} with
        filtered seqs
    """

    gc_bounds = fastq_tools.make_bounds(gc_bounds)
    length_bounds = fastq_tools.make_bounds(length_bounds)  # make tuple-like bounds

    passed_filtration_seqs = {}

    for read_name, (read_seq, read_quality) in seqs.items():
        if len(read_seq) == 0:  # dodge zero division error to a more understandable one
            raise ValueError("Cannnot work with sequence of length 0")

        has_passed_filters = (
            fastq_tools.check_if_in_bounds(
                fastq_tools.count_gc_content(read_seq), gc_bounds
            )  # is in gc bounds
            and fastq_tools.check_if_in_bounds(
                len(read_seq), length_bounds
            )  # in length bounds
            and fastq_tools.check_mean_quality(
                fastq_tools.count_mean_quality(read_quality), quality_threshold
            )  # quality greater than
        )

        if has_passed_filters:  # how can i avoid copy here?
            passed_filtration_seqs[read_name] = seqs[read_name]

    return passed_filtration_seqs


class FastQFilter:
    """
    A utility class for filtering FASTQ files based on quality, GC content, and read length.

    Args:
        path_to_input (str): Path to the input FASTQ file.
        path_to_output (str, optional): Path to the output filtered FASTQ file. If not provided,
            a unique output path will be generated based on the input file name.
        gc_bounds (Iterable[Number, Number] | Number, optional): GC content bounds for filtering.
            Can be a single value or a tuple representing the lower and upper bounds (default: (0, 100)).
        length_bounds (Iterable[Number, Number] | Number, optional): Read length bounds for filtering.
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

    Private Methods:
        _make_bounds(bounds) -> tuple: Converts bounds to a valid tuple format.
        _uniquify_path(path) -> str: Generates a unique output path.
        _is_in_bounds(value: Number, bounds: tuple) -> bool: Checks if a value is within bounds.
    """

    def __init__(
        self,
        path_to_input: str,
        path_to_output: str = None,
        gc_bounds: Iterable[Number, Number] | Number = (0, 100),
        length_bounds: Iterable[Number, Number] | Number = (0, 2**32),
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
        return f"FastQFilter(\n\tpath_to_input={self.path_to_input},\n\tpath_to_output={self.path_to_output}\n\tgc_bounds={self.gc_bounds},\n\tlength_bounds={self.length_bounds},\n\tquality_threshold={self.quality_threshold}\n)"


def run_ultimate_protein_tools(seqs: dict, command: str, **kwargs) -> dict:
    """
    Accepts command and runs it on input data with params

    Arguments:
    - seqs (str): Input in form of path, seq, seq list or seq dct
    - command (str): Valid command from command_dct
    - **kwargs to be passed to inner funcs

    Return:
    - output_dct (dict):
        dict where keys are number or name of seq and values are results of command run
    """

    output_dict = {}
    for seq_name, seq in seqs.items():
        if command in command_dct_prot:
            if command == "is_protein_valid":
                output_dict |= {seq_name: protein_tools.is_protein_valid(seq)}
            else:
                if protein_tools.is_protein_valid(seq):
                    output_dict |= {seq_name: command_dct_prot[command](seq, **kwargs)}
                else:
                    raise ValueError(f"Invalid protein, name/number: {seq_name}")

    if len(output_dict) == 1:
        return output_dict[list(output_dict.keys())[0]]
    return output_dict
