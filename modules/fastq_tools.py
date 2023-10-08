'''
Filtrate fastq reads
'''

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

