# TODO add module docstring
NUCL_COMP_DCT = {
    'RNA': {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
            'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'},
    'DNA': {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g'},
    'DNA_to_RNA': {'T': 'U', 'U': 'T', 't': 'u', 'u': 't'}
}


def check_seq_type(seq: str) -> str | None:
    '''
    Checks seq type (DNA, RNA or None). Returns str or None. Cases like ACG assumed to be DNA

    Arguments:
    - seq (str): given sequence

    Return:
    - seq_type (srt | None): Type of sequence as str or None
    '''

    dna = set(NUCL_COMP_DCT['DNA'])
    rna = set(NUCL_COMP_DCT['RNA'])
    seq = set(seq)

    seq_type = None

    if seq.issubset(rna) and not seq.issubset(dna): # in RNA and not in DNA
        seq_type = 'RNA'
    elif seq.issubset(dna): # in DNA and not in RNA (ACG (no T or U) cases considered DNA)
        seq_type ='DNA'
    return seq_type


def reverse(seq: str) -> str:
    '''
    Reverses given seq. 

    Arguments:
    - seq (str): given sequence

    Return:
    - srt: Reversed seq
    '''
    return seq[::-1]


def complement(seq: str, nucl_type: str) -> str:
    '''
    Complements given seq. Nucleic acid-type blind

    Arguments:
    - seq (str): given sequence

    Return:
    - srt: Complemented seq
    '''

    outseq = []
    for letter in seq:
        outseq.append(NUCL_COMP_DCT[nucl_type][letter])
    return ''.join(outseq)


def transcribe(seq: str) -> str:
    '''
    Transcribes given DNA to RNA or reverse transcribes RNA to DNA. Nucleic acid-type blind

    Arguments:
    - seq (str): given sequence

    Return:
    - str: Transcribed RNA seq
    '''

    outseq = []
    for letter in seq: # loop here to dodge O(n^2) double-replace case (can be sensitive on larger seq)
        if letter in NUCL_COMP_DCT['DNA_to_RNA']:
            letter = NUCL_COMP_DCT['DNA_to_RNA'][letter]
        outseq.append(letter)

    return ''.join(outseq)


def create_input_dict(*inp: str | list) -> dict:
    """
    Parses input seq or list of seqs and returns numerated dict of seqs.

    Arguments:
    - inp (str): Input seq or list of seqs

    Return:
    - parsed_dct (dict): numerated dict in format {0: 'seq'}
    """

    parsed_dct = {}
    if isinstance(inp, tuple):
        for i, seq in enumerate(inp):
            parsed_dct |= {i: seq}

    return parsed_dct