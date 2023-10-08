'''
Yet another "cool" bioinformatics tool to handle DNAs, RNAs, proteins and filter .fastq
'''

import modules.dna_rna_tools as dna_rna_tools
import modules.protein_tools as protein_tools
import modules.fastq_tools as fastq_tools


command_dict = {
    'check_seq_type': check_seq_type,
    'reverse': reverse,
    'complement': complement,
    'transcribe': transcribe,
}


def run_dna_rna_tools(seqs: dict, command: str) -> dict:
    '''
    Runs dna_rna_tools on given dict of seqs with given command

    Arguments:
    - seqs (dict): Input dict of format {seq_name: 'seq'}

    Return:
    - output_dict (dict | str | bool): 
        dict of results of operations {seq_name: 'result'} or one operation result 
        if one sequence given
    '''

    output_dict = {}

    for seq_name, seq in seqs.items():
        if command == 'check_seq_type': # user may want to check given seqs
            output_dict |= {seq_name: check_seq_type(seq)}
        else: # if other command
            nucl_type = check_seq_type(seq)
            if nucl_type is None:
                raise ValueError('Can only work with DNA or RNA sequence')

            if command == 'complement':
                output_dict |= {seq_name: complement(seq, nucl_type)}
            else:
                output_dict |= {seq_name: command_dict[command](seq)}

    if len(output_dict) == 1:
        return output_dict[list(output_dict.keys())[0]]
    return output_dict


