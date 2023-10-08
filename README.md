# YABT: yet another bioinformatic tool

This repository is the result of HW#5 on modules at Bioiformatics Institute.

Primary goal of this homework was to conclude recieved knowledge on functions and modules.

### Installation 

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


### Quickstart

1. General overview

YABT has three major functions available at `yabt` namespace which are:
- `yabt.run_dna_rna_tools` - to process DNA and RNA seqs
- `yabt.run_ultimate_protein_tools` - to process proteins
- `yabt.yabt.run_fastq_tools` - to filter `fastq`-reads

2. Input

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

3. `yabt.run_dna_rna_tools` overview

This module can perform the following operations on sequences:

- `'check_seq_type'`

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

- `'reverse'`

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

- `'complement'`

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

- `'transcribe'`

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

