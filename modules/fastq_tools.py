def count_gc_content(seq: str) -> float:
    '''
    Counts GC-content in a given sequence

    Arguments:
    - seq (str): given sequence

    Return:
    - gc_content (float): GC-ratio of a given seq as a percentage
    '''

    gc_count = 0
    for letter in seq.upper(): # used loop instead of .count method to iterate over seq only once (O(n))
        if letter == 'G' or letter == 'C':
            gc_count += 1
    if len(seq) == 0:
        raise ValueError('Cannnot work with sequence of length 0') # dodge zero division error to a more understandable one
    gc_content = gc_count / len(seq) * 100 # len() is believed to be O(1)
    return gc_content


def make_bounds(bounds: int | float | list | tuple) -> tuple:
    '''
    Checks input bounds to be lists or tuple or converts bounds from int and float to format (0, int) or (0, float)

    Arguments:
    - bounds (int | float | list | tuple): given bounds

    Return:
    - bounds (tuple): bounds in format of tuple
    '''
    
    if type(bounds) == int or type(bounds) == float:
        return (0, bounds)
    elif type(bounds) == list or type(bounds) == tuple:
        return tuple(bounds) # make return more expexted for further maintaining if list passed
    else:
        raise TypeError(f'Cannot work with {type(bounds).__name__} type')
    

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
