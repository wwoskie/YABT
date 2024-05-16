import os

from dataclasses import dataclass


@dataclass
class FastaRecord:
    """
    Represents a FASTA record containing sequence information.

    Attributes:
        id (str): Identifier for the sequence.
        description (str): Description or additional information about the sequence.
        seq (str): The actual sequence data.
    """

    id: str
    description: str
    seq: str

    def __repr__(self):
        return (
            f"FastaRecord(\n\tid='{self.id}'\n\t"
            + f"description='{self.description}'"
            + f"\n\tseq={self.seq}\n)"
        )


class OpenFasta:
    """
    Context manager for reading FASTA files.

    Args:
        path_to_fasta (str): Path to the FASTA file.

    Methods:
        read_record(): Reads the next FASTA record and returns a FastaRecord object.
        read_records(): Reads all FASTA records and returns a list of FastaRecord objects.
    """

    def __init__(self, path_to_fasta):
        self.path_to_fasta = path_to_fasta
        self.handler = None
        self.current_info = None

    def __enter__(self):
        self.handler = open(self.path_to_fasta)
        self.current_info = self.handler.readline().strip()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.handler:
            self.handler.close()

    def __iter__(self):
        return self

    def __next__(self) -> FastaRecord:
        next_info = None
        seq_list = []

        for line in self.handler:
            line = line.strip()
            if line.startswith(">"):
                next_info = line
                break
            else:
                seq_list.append(line)

        if not next_info:
            next_info = self.current_info

        splitted_line = self.current_info[1:].split(" ")
        seq_id = splitted_line.pop(0)
        seq_description = " ".join(splitted_line)

        if seq_list == []:
            raise StopIteration

        seq = "".join(seq_list)
        self.current_info = next_info

        return FastaRecord(id=seq_id, description=seq_description, seq=seq)

    def read_record(self) -> FastaRecord:
        """Reads the next FASTA record and returns a FastaRecord object."""

        return self.__next__()

    def read_records(self) -> list[FastaRecord]:
        """Reads all FASTA records and returns a list of FastaRecord objects."""

        records = []

        for record in self:
            records.append(record)

        return records
