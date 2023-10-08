'''
Manupulate protein sequences.
'''

RNA_AA_TABLE = {
    'F': ['UUU', 'UUC'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'Y': ['UAU', 'UAC'],
    '*': ['UAA', 'UAG', 'UGA', 'uaa', 'uag', 'uga'],
    'C': ['UGU', 'UGC'],
    'W': ['UGG'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'H': ['CAU', 'CAC'],
    'Q': ['CAA', 'CAG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'I': ['AUU', 'AUC', 'AUA'],
    'M': ['AUG'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'N': ['AAU', 'AAC'],
    'K': ['AAA', 'AAG'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'D': ['GAU', 'GAC'],
    'E': ['GAA', 'GAG'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'f': ['uuu', 'uuc'],
    'l': ['uua', 'uug', 'cuu', 'cuc', 'cua', 'cug'],
    's': ['ucu', 'ucc', 'uca', 'ucg', 'agu', 'agc'],
    'y': ['uau', 'uac'],
    'c': ['ugu', 'ugc'],
    'w': ['ugg'],
    'p': ['ccu', 'ccc', 'cca', 'ccg'],
    'h': ['cau', 'cac'],
    'q': ['caa', 'cag'],
    'r': ['cgu', 'cgc', 'cga', 'cgg', 'aga', 'agg'],
    'i': ['auu', 'auc', 'aua'],
    'm': ['aug'],
    't': ['acu', 'acc', 'aca', 'acg'],
    'n': ['aau', 'aac'],
    'k': ['aaa', 'aag'],
    'v': ['guu', 'guc', 'gua', 'gug'],
    'a': ['gcu', 'gcc', 'gca', 'gcg'],
    'd': ['gau', 'gac'],
    'e': ['gaa', 'gag'],
    'g': ['ggu', 'ggc', 'gga', 'ggg']
}

RNA_CODON_TABLE = { 
  'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'UCU': 'S',
  'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'UAU': 'Y', 'UAC': 'Y',
  'UAA': '*', 'UAG': '*', 'UGU': 'C', 'UGC': 'C', 'UGA': '*',
  'UGG': 'W', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
  'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAU': 'H',
  'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R', 'CGC': 'R',
  'CGA': 'R', 'CGG': 'R', 'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
  'AUG': 'M', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
  'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGU': 'S',
  'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V', 'GUC': 'V',
  'GUA': 'V', 'GUG': 'V', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A',
  'GCG': 'A', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
  'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'uuu': 'f',
  'uuc': 'f', 'uua': 'l', 'uug': 'l', 'ucu': 's', 'ucc': 's',
  'uca': 's', 'ucg': 's', 'uau': 'y', 'uac': 'y', 'uaa': '*',
  'uag': '*', 'ugu': 'c', 'ugc': 'c', 'uga': '*', 'ugg': 'w',
  'cuu': 'l', 'cuc': 'l', 'cua': 'l', 'cug': 'l', 'ccu': 'p',
  'ccc': 'p', 'cca': 'p', 'ccg': 'p', 'cau': 'h', 'cac': 'h',
  'caa': 'q', 'cag': 'q', 'cgu': 'r', 'cgc': 'r', 'cga': 'r',
  'cgg': 'r', 'auu': 'i', 'auc': 'i', 'aua': 'i', 'aug': 'm',
  'acu': 't', 'acc': 't', 'aca': 't', 'acg': 't', 'aau': 'n',
  'aac': 'n', 'aaa': 'k', 'aag': 'k', 'agu': 's', 'agc': 's',
  'aga': 'r', 'agg': 'r', 'guu': 'v', 'guc': 'v', 'gua': 'v',
  'gug': 'v', 'gcu': 'a', 'gcc': 'a', 'gca': 'a', 'gcg': 'a',
  'gau': 'd', 'gac': 'd', 'gaa': 'e', 'gag': 'e', 'ggu': 'g',
  'ggc': 'g', 'gga': 'g', 'ggg': 'g'
}


def read_seq_from_fasta(path_to_seq: str,
                        use_full_name: bool = False,
                        **_) -> dict:
    """
    Reads sequences from fasta file and returns dictionary.

    Arguments:
    - path_to_seq (str): path to file

    Return:
    - dict: dict of sequences names as keys and sequences themselves as values {'seq_name': 'sequence',}
    """

    with open(path_to_seq) as f:
        out_dct = {}
        name = None # None is set to skip first ''.join(current_seqs) 
        for line in f:
            line = line.strip()
            if line.startswith('>'):  # check for first line in seq
                if not name is None:
                    out_dct[name] = ''.join(current_seqs) # join current_seqs to str if not first seq_name
                if use_full_name:  # check if user set full name in fasta
                    name = line[1:]  # take whole fasta properties (e.g. if names not unique)
                else:
                    name = line[1:].split()[0]
                current_seqs = [] # create list to append read lines to
            else:
                current_seqs.append(line)  # get value from dict (return '' if empty) and append str

        out_dct[name] = ''.join(current_seqs) # collects last seqs to dict
            
    return out_dct


def get_sites_lengths(sites: list) -> dict:
    """
    Takes sites list and calculates their lengths. Used inside find_sites func

    Arguments:
    - sites (list): list of sites (str)

    Return:
    - dict: dict of sites length {'site': 'length',}
    """

    sites_length_dct = {}
    for site in sites:
        sites_length_dct[site] = len(site)
    return sites_length_dct


def invert_dct(dct: dict) -> dict:
    """
    Inverts a dict. Used inside find_sites func

    Arguments:
    - dct (dict): dict to be inverted

    Return:
    - dict: inverted dict
    """

    inv_dct = {}
    for k, v in dct.items():
        inv_dct[v] = inv_dct.get(v, []) + [k]  # get value from dict (return [] if empty) and append key
    return inv_dct
