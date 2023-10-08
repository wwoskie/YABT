# TODO add module docstring
NUCL_COMP_DCT = {'RNA': {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
                         'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'},
                 'DNA': {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                         'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}}


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
    Reverses given seq

    Arguments:
    - seq (str): given sequence

    Return:
    - srt: Reversed seq
    '''
    return seq[::-1]