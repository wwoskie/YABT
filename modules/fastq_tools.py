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
