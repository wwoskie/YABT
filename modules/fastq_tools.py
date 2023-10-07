def count_gc_content(seq: str):
    '''
    Counts GC-content in a given sequence

    Arguments:
    - seq (str): given sequence

    Return:
    - float: GC-ratio of a given seq as a percentage
    '''

    gc_count = 0
    for letter in seq.upper(): # used loop instead of .count method to iterate over seq only once (O(n))
        if letter == 'G' or letter == 'C':
            gc_count += 1
    gc_content = gc_count / len(seq) * 100 # len() is believed to be O(1)
    return gc_content
