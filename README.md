# YABT: yet another bioinformatic tool

This repository is the result of HW#5 on modules at Bioiformatics Institute.

Primary goal of this homework was to conclude recieved knowledge on functions and modules.

## Installation 

This tool can be installed via the following ways:
- `git clone https://github.com/wwoskie/YABT.git` to either local path of your script (`.` directory) and import it directly with:

```python 
import yet_another_bioinformatic_tool as yabt
```

- or you can choose to specify location during imports with:
```python
import sys
sys.path.append('/path/to/yet_another_bioinformatic_tool.py')
import yet_another_bioinformatic_tool as yabt
```
One can use any custom name for import but we suggest using `yabt` as module shortname.

But don't limit yourself to only these two ways, as you can find more on modules as a whole [here](https://docs.python.org/3/tutorial/modules.html)


## Quickstart


### 1. General overview

YABT has three major functions available at `yabt` namespace which are:
- `yabt.run_dna_rna_tools` - to process DNA and RNA seqs
- `yabt.run_ultimate_protein_tools` - to process proteins
- `yabt.yabt.run_fastq_tools` - to filter `fastq`-reads


### 2. Input

For these functions dictionary is the only supported input type. However, `yabt.yabt.run_fastq_tools` reqires slightly different input dictionary of format `{'seq_name': ('nucl_seq', 'quality_for_seq')}` while other two use `{'seq_name': 'seq'}` type of dictionary.

To create this type of input dict for `yabt.run_dna_rna_tools` and `yabt.run_ultimate_protein_tools` two functions can be used:

- `yabt.create_input_dict`:

Parses input seq or list of seqs and returns numerated dict of seqs.

    Arguments:
    - inp (str): Input seq or list of seqs

    Return:
    - parsed_dct (dict): numerated dict in format {0: 'seq'}

Examples:
```python
yabt.create_input_dict('AUGC', 'acga')
> {0: 'AUGC', 1: 'acga'}
```

- `yabt.read_seq_from_fasta`:

Reads sequences from fasta file and returns dictionary.

    Arguments:
    - path_to_seq (str): path to file

    Return:
    - dict: 
        dict of sequences names as keys and sequences themselves as values {'seq_name': 'sequence'}

Examples:
```python
fasta_seqs = yabt.read_seq_from_fasta('testdata/testdata.fasta')
fasta_seqs_full = yabt.read_seq_from_fasta('testdata/testdata.fasta', use_full_name=True)
print(fasta_seqs)
> {'crab_anapl': 'MDITIHNPLIRRPLFSWLAPSRIFDQIFGEHLQESELLPASPSLSPFLMRSPIFRMPSWLETGLSEMRLEKDKFSVNLDVKHFSPEELKVKVLGDMVEIHGKHEERQDEHGFIAREFNRKYRIPADVDPLTITSSLSLDGVLTVSAPRKQSDVPERSIPITREEKPAIAGAQRK', 'crab_bovin': 'MDIAIHHPWIRRPFFPFHSPSRLFDQFFGEHLLESDLFPASTSLSPFYLRPPSFLRAPSWIDTGLSEMRLEKDRFSVNLDVKHFSPEELKVKVLGDVIEVHGKHEERQDEHGFISREFHRKYRIPADVDPLAITSSLSSDGVLTVNGPRKQASGPERTIPITREEKPAVTAAPKK', 'crab_chick': 'MDITIHNPLVRRPLFSWLTPSRIFDQIFGEHLQESELLPTSPSLSPFLMRSPFFRMPSWLETGLSEMRLEKDKFSVNLDVKHFSPEELKVKVLGDMIEIHGKHEERQDEHGFIAREFSRKYRIPADVDPLTITSSLSLDGVLTVSAPRKQSDVPERSIPITREEKPAIAGSQRK', 'crab_human': 'MDIAIHHPWIRRPFFPFHSPSRLFDQFFGEHLLESDLFPTSTSLSPFYLRPPSFLRAPSWFDTGLSEMRLEKDRFSVNLDVKHFSPEELKVKVLGDVIEVHGKHEERQDEHGFISREFHRKYRIPADVDPLTITSSLSSDGVLTVNGPRKQVSGPERTIPITREEKPAVTAAPKK', 'crab_mesau': 'MDIAIHHPWIRRPFFPFHSPSRLFDQFFGEHLLESDLFSTATSLSPFYLRPPSFLRAPSWIDTGLSEMRMEKDRFSVNLDVKHFSPEELKVKVLGDVVEVHGKHEERQDEHGFISREFHRKYRIPADVDPLTITSSLSSDGVLTVNGPRKQASGPERTIPITREEKPAVTAAPKK', 'crab_mouse': 'MDIAIHHPWIRRPFFPFHSPSRLFDQFFGEHLLESDLFSTATSLSPFYLRPPSFLRAPSWIDTGLSEMRLEKDRFSVNLDVKHFSPEELKVKVLGDVIEVHGKHEERQDEHGFISREFHRKYRIPADVDPLAITSSLSSDGVLTVNGPRKQVSGPERTIPITREEKPAVAAAPKK', 'crab_rabit': 'MDIAIHHPWIRRPFFPFHSPSRLFDQFFGEHLLESDLFPTSTSLSPFYLRPPSFLRAPSWIDTGLSEMRLEKDRFSVNLDVKHFSPEELKVKVLGDVIEVHGKHEERQDEHGFISREFHRKYRIPADVDPLTITSSLSSDGVLTVNGPRKQAPGPERTIPITREEKPAVTAAPKK', 'crab_rat': 'MDIAIHHPWIRRPFFPFHSPSRLFDQFFGEHLLESDLFSTATSLSPFYLRPPSFLRAPSWIDTGLSEMRMEKDRFSVNLDVKHFSPEELKVKVLGDVIEVHGKHEERQDEHGFISREFHRKYRIPADVDPLTITSSLSSDGVLTVNGPRKQASGPERTIPITREEKPAVTAAPKK', 'crab_squac': 'MDIAIQHPWLRRPLFPSSIFPSRIFDQNFGEHFDPDLFPSFSSMLSPFYWRMGAPMARMPSWAQTGLSELRLDKDKFAIHLDVKHFTPEELRVKILGDFIEVQAQHEERQDEHGYVSREFHRKYKVPAGVDPLVITCSLSADGVLTITGPRKVADVPERSVPISRDEKPAVAGPQQK'}
```


### 3. `yabt.run_dna_rna_tools` overview

Output of operation of the module will be unpacked if only one seq is given. This module can perform the following operations on nucleic sequences:

#### - `'check_seq_type'`

Checks seq type (DNA, RNA or None). Returns str or None. Cases like ACG assumed to be DNA

    Arguments:
    - seq (str): given sequence

    Return:
    - seq_type (srt | None): Type of sequence as str or None

Can be useful if user is not sure in quality of processed data, Returns sequence type as `str` `'DNA'` or `'RNA'` or `None` if seq is corrupt.

Examples:
```python
yabt.run_dna_rna_tools(yabt.create_input_dict('AUGC'), 'check_seq_type')
> 'RNA'
yabt.run_dna_rna_tools(yabt.create_input_dict('ATGC', 'cuga', 'kek'), 'check_seq_type')
> {0: 'DNA', 1: 'RNA', 2: None}
yabt.run_dna_rna_tools({'my_cool_DNA': 'TGCa', 'my_cool_RNA': 'uaaugauga', 'lol': 'kek'}, 'check_seq_type')
> {'my_cool_DNA': 'DNA', 'my_cool_RNA': 'RNA', 'lol': None}
print(yabt.run_dna_rna_tools(yabt.create_input_dict('AutGC'), 'check_seq_type'))
print(yabt.run_dna_rna_tools(yabt.create_input_dict('I love Python'), 'check_seq_type'))
> None
> None
```

#### - `'reverse'`

Reverses given seq. 

    Arguments:
    - seq (str): given sequence

    Return:
    - srt: Reversed seq
    
Simple reverse operation that, however, checks if seq is valid and throws ValueError if not

Examples:
```python
yabt.run_dna_rna_tools(yabt.create_input_dict('ATGC'), 'reverse')
> 'CGTA'
yabt.run_dna_rna_tools(yabt.create_input_dict('ATGC', 'ctga'), 'reverse')
> {0: 'CGTA', 1: 'agtc'}
yabt.run_dna_rna_tools(yabt.create_input_dict('lol'), 'reverse')
> ... ValueError: Can only work with DNA or RNA sequence
```

#### - `'complement'`

Complements given seq. Nucleic acid-type blind

    Arguments:
    - seq (str): given sequence
    - nucl_type (str): type of nucleic acid

    Return:
    - srt: Complemented seq

Complements the sequence. Is nucleic acid-type blind (though requires nucleic acid type as input), so can potentionally be used for RNA secondary structures prediction (though has limited pair options by default). Complement dictionary can be modified though by experienced user.

Examples:
```python
yabt.run_dna_rna_tools(yabt.create_input_dict('ATGC', 'cuga'), 'complement')
> {0: 'TACG', 1: 'gacu'}
```

#### - `'transcribe'`

Transcribes given DNA to RNA or reverse transcribes RNA to DNA. Nucleic acid-type blind

    Arguments:
    - seq (str): given sequence

    Return:
    - str: Transcribed RNA seq
    '''

This function can perform not only transcription, but reverse-transcription too as it is nucleic acid-type blind.

Examples:
```python
yabt.run_dna_rna_tools(yabt.create_input_dict('ATGC', 'cuga'), 'transcribe')
> {0: 'AUGC', 1: 'ctga'}
yabt.run_dna_rna_tools(yabt.create_input_dict('ATtAcGC'), 'transcribe')
> 'AUuAcGC'
```

All source functions are located at `./modules/dna_rna_tools.py` and can be modified by user for their needs. Also note that this module supports masking and keeps lower/uppercase during operations.

If you need any internal individual function, you can access it by importing it directly 
```python
from yabt.dna_rna_tools import <function_name> as <your_function_name>
```
This also applies for all further modules


### 4. `yabt.run_protein_tools` overview

This module can perform the following operations on protein sequences:

#### - `'find_sites'`

Finds indexes of given sites.

    Arguments:
    - seq (str): seq to be checked
    - sites (list): sites to be found in form of a list
    - is_one_based (bool): whether result should be 0- (False) or 1-indexed (True). Default False

    Return:
    - dict: dictionary of sites as keys and lists of indexes for the site where it's been found

Finds all indexes of given sites in a seq. Note `is_one_based` parameter. Please, don't function output with this parameter further in python (0-index-based language) to avoid any mess. Though it might be useful for human or R input.

Examples:
```python
yabt.run_dna_rna_tools(yabt.create_input_dict('ATGC', 'cuga'), 'transcribe')
> {0: 'AUGC', 1: 'ctga'}
fasta_seqs = yabt.read_seq_from_fasta('/content/testdata.fasta')
yabt.run_ultimate_protein_tools(fasta_seqs, 'find_sites', sites = ['M', 'MD', 'RPLF'])
> {'crab_anapl': {'M': [0, 48, 55, 66, 95], 'MD': [0], 'RPLF': [11]},
 'crab_bovin': {'M': [0, 67], 'MD': [0]},
 'crab_chick': {'M': [0, 48, 55, 66, 95], 'MD': [0], 'RPLF': [11]},
 'crab_human': {'M': [0, 67], 'MD': [0]},
 'crab_mesau': {'M': [0, 67, 69], 'MD': [0]},
 'crab_mouse': {'M': [0, 67], 'MD': [0]},
 'crab_rabit': {'M': [0, 67], 'MD': [0]},
 'crab_rat': {'M': [0, 67, 69], 'MD': [0]},
 'crab_squac': {'M': [0, 43, 51, 55, 58], 'MD': [0], 'RPLF': [11]}}
```

#### - `'get_protein_rnas_number'`

Get number of all possible RNA's for a given protein.

    Arguments:
    - seq (str): seq to be checked

    Return:
    - rnas_num (int): number of possible RNA's for seq

Examples:
```python
yabt.run_ultimate_protein_tools({'my_seq': 'AAAAAAAAA', 'my_seq2': 'NNnnNN'}, 'get_protein_rnas_number')
> {'my_seq': 262144, 'my_seq2': 64}
```

#### - `'is_protein_valid'`

Checks if protein is valid.

    Arguments:
    - seq (str): seq to be checked

    Return:
    - bool, the result of the check

Can be used for data checking if user not sure in input data

Examples:
```python
yabt.run_ultimate_protein_tools({'my_seq': 'KEK', 'my_seq2': 'lol'}, 'is_protein_valid')
> {'my_seq': True, 'my_seq2': False}
```

#### - `'get_length_of_protein'`

Calculates the length of a protein.

    Arguments:
    - seq (str): sequence to calculate the length

    Return:
    - int: sequence length

Also checks if protein is valid and throws ValueError if not

Examples:
```python
yabt.run_ultimate_protein_tools({'kek_seq': 'KEK', 'lol_seq': 'LLLLL'}, 'get_length_of_protein')
> {'kek_seq': 3, 'lol_seq': 5}
yabt.run_ultimate_protein_tools({'kek_seq': 'KEK', 'lol_seq': 'lol'}, 'get_length_of_protein')
> ... ValueError: Invalid protein, name/number: lol_seq
```

#### -`'count_aa'`

Counts the number of given or all amino acids in a protein sequence.

    Arguments:
    - seq (str): sequence to count amino acids
    - aminoacids (str): which amino acids to count in sequence

    Return:
    - dict: a dictionary with amino acids and its count

Examples:
```python
yabt.run_ultimate_protein_tools(fasta_seqs, 'count_aa')['crab_rat']
> {'M': 3,
 'D': 11,
 'I': 10,
 'A': 8,
 'H': 9,
 'P': 16,
 'W': 2,
 'R': 14,
 'F': 13,
 'S': 17,
 'L': 14,
 'Q': 3,
 'G': 8,
 'E': 14,
 'T': 9,
 'Y': 2,
 'K': 10,
 'V': 10,
 'N': 2}
yabt.run_ultimate_protein_tools(fasta_seqs, 'get_fracture_of_aa', aminoacids_to_count='ML')['crab_rat']
> {'M': 0.0171, 'L': 0.08}
yabt.run_ultimate_protein_tools(fasta_seqs, 'count_aa', aminoacids_to_count='ARS')['crab_rat']
{'A': 8, 'R': 14, 'S': 17}
```

#### - `'get_fracture_of_aa'`

Calculates the fracture or percentage of amino acids in a protein sequence.

    Arguments:
    - seq (str): sequence in which you need to calculate the fracture of amino acids
    - show_as_percentage (bool): change it to True, if you want to get results with percentages
    - aminoacids (str): the fracture of which amino acids to count in the sequence

    Return:
    - dict: a dictionary with amino acids and its fracture or percentage

Examples:
```python
yabt.run_ultimate_protein_tools(fasta_seqs, 'get_fracture_of_aa')['crab_rat']
> {'M': 0.0171,
 'D': 0.0629,
 'I': 0.0571,
 'A': 0.0457,
 'H': 0.0514,
 'P': 0.0914,
 'W': 0.0114,
 'R': 0.08,
 'F': 0.0743,
 'S': 0.0971,
 'L': 0.08,
 'Q': 0.0171,
 'G': 0.0457,
 'E': 0.08,
 'T': 0.0514,
 'Y': 0.0114,
 'K': 0.0571,
 'V': 0.0571,
 'N': 0.0114}}
yabt.run_ultimate_protein_tools(fasta_seqs, 'get_fracture_of_aa', aminoacids_to_count='ML')['crab_rat']
> {'M': 0.0171, 'L': 0.08}
```