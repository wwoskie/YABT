'''
Yet another "cool" bioinformatics tool to handle DNAs, RNAs, proteins and filter .fastq
'''

import modules.dna_rna_tools as dna_rna_tools
import modules.protein_tools as protein_tools
import modules.fastq_tools as fastq_tools


command_dict_nucl = {
    'check_seq_type': dna_rna_tools.check_seq_type,
    'reverse': dna_rna_tools.reverse,
    'complement': dna_rna_tools.complement,
    'transcribe': dna_rna_tools.transcribe,
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
            output_dict |= {seq_name: dna_rna_tools.check_seq_type(seq)}
        else: # if other command
            nucl_type = dna_rna_tools.check_seq_type(seq)
            if nucl_type is None:
                raise ValueError('Can only work with DNA or RNA sequence')

            if command == 'complement':
                output_dict |= {seq_name: dna_rna_tools.complement(seq, nucl_type)}
            else:
                output_dict |= {seq_name: command_dict[command](seq)}

    if len(output_dict) == 1:
        return output_dict[list(output_dict.keys())[0]]
    return output_dict


def run_fastq_tools(seqs: dict, # how can i make native type hint here?
                    gc_bounds: tuple | list| int | float = (0, 100),
                    length_bounds: tuple | list| int | float = (0, 2**32),
                    quality_threshold: int | float = 0) -> dict:
    '''
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
    '''

    gc_bounds = fastq_tools.make_bounds(gc_bounds)
    length_bounds = fastq_tools.make_bounds(length_bounds) # make tuple-like bounds

    passed_filtration_seqs = {}

    for read_name, (read_seq, read_quality) in seqs.items(): # TODO add is_dna check
        if len(read_seq) == 0: # dodge zero division error to a more understandable one
            raise ValueError('Cannnot work with sequence of length 0')

        has_passed_filters = (
            fastq_tools.check_if_in_bounds(fastq_tools.count_gc_content(read_seq),
                               gc_bounds) # is in gc bounds
            and fastq_tools.check_if_in_bounds(len(read_seq),
                                   length_bounds) # in length bounds
            and fastq_tools.check_mean_quality(fastq_tools.count_mean_quality(read_quality),
                                   quality_threshold) # quality greater than
            )

        if has_passed_filters: # how can i avoid copy here?
            passed_filtration_seqs[read_name] = seqs[read_name]

    return passed_filtration_seqs