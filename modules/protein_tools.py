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
    - dict: 
        dict of sequences names as keys and sequences themselves as values {'seq_name': 'sequence',}
    """

    with open(path_to_seq, 'r') as f:
        out_dct = {}
        name = None # None is set to skip first ''.join(current_seqs)
        for line in f:
            line = line.strip()
            if line.startswith('>'):  # check for first line in seq
                if not name is None: # join current_seqs to str if not first seq_name
                    out_dct[name] = ''.join(current_seqs)
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
    for k, v in dct.items(): # get value from dict (return [] if empty) and append key
        inv_dct[v] = inv_dct.get(v, []) + [k]
    return inv_dct


def is_protein_valid(seq: str) -> bool:
    """
    Checks if protein is valid.

    Arguments:
    - seq (str): seq to be checked

    Return:
    - bool, the result of the check
    """

    if set(seq).issubset(RNA_AA_TABLE):
        return True
    return False


def find_sites(seq: str,
               sites: list,
               is_one_based: bool = False
               ) -> dict:
    """
    Finds indexes of given sites.

    Arguments:
    - seq (str): seq to be checked
    - sites (list): sites to be found in form of a list
    - is_one_based (bool): whether result should be 0- (False) or 1-indexed (True). Default False

    Return:
    - dict: dictionary of sites as keys and lists of indexes for the site where it's been found
    """

    window_sizes = invert_dct(get_sites_lengths(sites))
    # get lengths of all sites and stick them together to avoid
    # passing through seq multiple times if possible
    found_sites = {}
    # perform iteration for all given lengths of sites
    # pls don't be scared of nested loops 2 of them are supposed to be relatively small
    # if user is not intended to search huge amount of sites
    for window_size, sites_of_window_size in window_sizes.items():
        for i in range(len(seq) - window_size + 1):
            # iterate through seq with step one and consider window
            # of site length each step
            scatter = seq[i:i + window_size]
            # get fragment of sequence with length of window i.e. scatter
            for site in sites_of_window_size:
                if scatter == site:  # check if scatter is site
                    found_sites[site] = (
                            found_sites.get(site, [])  # get
                            + [i + is_one_based]
                    )  # append index to list in dict
    return found_sites


def get_protein_rnas_number(seq: int, **_) -> int:
    """
    Get number of all possible RNA's for a given protein.

    Arguments:
    - seq (str): seq to be checked

    Return:
    - rnas_num (int): number of possible RNA's for seq
    """

    rnas_num = 1
    for amino_acid in seq:
        rnas_num *= len(RNA_AA_TABLE[amino_acid])
    return rnas_num


def get_length_of_protein(seq: str) -> int:
    """
    Calculates the length of a protein.

    Arguments:
    - seq (str): sequence to calculate the length

    Return:
    - int: sequence length
    """

    return len(seq)


def count_aa(seq: str, aminoacids_to_count: str = None) -> dict:
    """
    Counts the number of given or all amino acids in a protein sequence.

    Arguments:
    - seq (str): sequence to count amino acids
    - aminoacids (str): which amino acids to count in sequence

    Return:
    - dict: a dictionary with amino acids and its count
    """

    aa_dict_count = {}
    if (aminoacids_to_count is None) or (aminoacids_to_count == ''):
        '''
        I added an additional condition for user-friendly experience.
        E.g., we can want to find specific aminoacid, look on result and then look on all aminoacids.
        Without this condition we have to delete keyword argument, but with it we can only make it empty.
        '''
        aminoacids_to_count = set(seq)
    for aa in seq:
        if aa in aminoacids_to_count:
            aa_dict_count[aa] = aa_dict_count.get(aa, 0) + 1 # return 1 if aa met first time, add 1 if not
    return aa_dict_count

def get_fracture_of_aa(seq: str, show_as_percentage: bool = False, aminoacids: str = None) -> dict:
    """
    Calculates the fracture or percentage of amino acids in a protein sequence.

    Arguments:
    - seq (str): sequence in which you need to calculate the fracture of amino acids
    - show_as_percentage (bool): change it to True, if you want to get results with percentages
    - aminoacids (str): the fracture of which amino acids to count in the sequence

    Return:
    - dict: a dictionary with amino acids and its fracture or percentage
    """

    if show_as_percentage:
        mult = 100
        round_var = 2
    else:
        mult = 1
        round_var = 4
    aa_dict_count = count_aa(seq, aminoacids=aminoacids)
    aa_dict_percent = {}
    len_of_protein = get_length_of_protein(seq)
    for aa, count in aa_dict_count.items():
        aa_dict_percent[aa] = round(count / len_of_protein * mult, round_var)
    return aa_dict_percent

