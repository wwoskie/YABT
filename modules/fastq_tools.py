def count_gc_content(seq: str) -> float:
    '''
    Counts GC-content in a given sequence

    Arguments:
    - seq (str): given sequence

    Return:
    - gc_content (float): GC-ratio of a given seq as a percentage
    '''

    gc_count = 0
    for letter in seq.upper(): # loop to avoid double-checking
        if letter in 'GC':
            gc_count += 1
    gc_content = gc_count / len(seq) * 100 # len() is believed to be O(1)
    return gc_content


def make_bounds(bounds: int | float | list | tuple) -> tuple:
    '''
    Checks input bounds to be lists or tuple. 
    Or converts bounds from int and float to format (0, int) or (0, float)

    Arguments:
    - bounds (int | float | list | tuple): given bounds

    Return:
    - bounds (tuple): bounds in format of tuple
    '''

    if not isinstance(bounds, (int, float, tuple, list)):
        raise TypeError(f'Cannot work with {type(bounds).__name__} type')

    if isinstance(bounds, (int, float)):
        bounds = (0, bounds)
    else:
        bounds = tuple(bounds) # make return more expected for further maintaining if list passed

    return bounds


def check_if_in_bounds(value: int | float, bounds: tuple) -> bool:
    '''
    Checks if value is in bounds (start and end included!) and returns True/False

    Arguments:
    - value (int | float): given value
    - bounds (tuple): bounds

    Return:
    - (bool): if value in bounds
    '''

    return bounds[0] <= value <= bounds[1]


def count_mean_quality(quality_seq: str) -> float:
    '''
    Counts mean read quality for a given quality str

    Arguments:
    - quality_seq (str): given quality str in ASCII format

    Return:
    - mean_quality (float): mean quality of read
    '''

    q_score_lst = []

    for single_letter_quality in quality_seq:
        q_score = ord(single_letter_quality) - 33
        q_score_lst.append(q_score)
    mean_quality = sum(q_score_lst) / len(q_score_lst)
    return mean_quality

def check_mean_quality(mean_quality: float, quality_threshold: int | float) -> bool:
    '''
    Counts if mean read quality is greater or equal to quality_threshold

    Arguments:
    - mean_quality (float): mean quality
    - quality_threshold (int | float): quality threshold to check against

    Return:
    - (bool): if mean_quality is greater or equal to quality_threshold
    '''

    return mean_quality >= quality_threshold


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

    gc_bounds = make_bounds(gc_bounds)
    length_bounds = make_bounds(length_bounds) # make tuple-like bounds

    passed_filtration_seqs = {}

    for read_name, (read_seq, read_quality) in seqs.items(): # TODO add is_dna check
        if len(read_seq) == 0: # dodge zero division error to a more understandable one
            raise ValueError('Cannnot work with sequence of length 0')

        has_passed_filters = (
            check_if_in_bounds(count_gc_content(read_seq),
                               gc_bounds) # is in gc bounds
            and check_if_in_bounds(len(read_seq),
                                   length_bounds) # in length bounds
            and check_mean_quality(count_mean_quality(read_quality),
                                   quality_threshold) # quality greater than
            )

        if has_passed_filters: # how can i avoid copy here?
            passed_filtration_seqs[read_name] = seqs[read_name]

    return passed_filtration_seqs
