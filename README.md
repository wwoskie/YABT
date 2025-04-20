# YABT: yet another bioinformatic tool

This repository is the result of HW#5 on modules at Bioiformatics Institute.

Primary goal of this homework was to conclude recieved knowledge on functions and modules.

## Installation 

This tool can be installed via the following ways:
- `git clone https://github.com/wwoskie/YABT.git` to either local path of your script (`.` directory) and import it directly with:

- or you can choose to specify location with:
```python
import sys
sys.path.append('/path/to/yet_another_bioinformatic_tool.py')
```
Import package classes via:
```python 
from yet_another_bioinformatic_tool import (
    AminoAcidSequence,
    DNASequence,
    FastQFiltrator,
    NucleicAcidSequence,
)
```

One can use any custom name for import but we suggest using `yabt` as module shortname.

But don't limit yourself to only these two ways, as you can find more on modules as a whole [here](https://docs.python.org/3/tutorial/modules.html)

## Requirements

Requirements can be found at `requirements.txt`

## Usage

Usage examples can be found at `examples.ipynb`


Alternatively, fastq filtration functionality can be accessed via CLI interface:

`python filter_fastq.py -i example_data/example_fastq.fastq -o output.fastq`

```bash
>>> python filter_fastq.py --help                                                
usage: filter_fastq.py [-h] [--input_path INPUT_PATH] [--output_path OUTPUT_PATH]
                       [--gc_lower_bound GC_LOWER_BOUND] [--gc_upper_bound GC_UPPER_BOUND]
                       [--length_lower_bound LENGTH_LOWER_BOUND]
                       [--length_upper_bound LENGTH_UPPER_BOUND]
                       [--quality_threshold QUALITY_THRESHOLD] [--logs_dir LOGS_DIR] [--rewrite]

Filter fastq files

options:
  -h, --help            show this help message and exit
  --input_path, -i INPUT_PATH
                        Path to input (default: None)
  --output_path, -o OUTPUT_PATH
                        Path to output (default: None)
  --gc_lower_bound, -gcl GC_LOWER_BOUND
                        GC content lower bound (default: 0)
  --gc_upper_bound, -gcu GC_UPPER_BOUND
                        GC content upper bound (default: 1)
  --length_lower_bound, -ll LENGTH_LOWER_BOUND
                        Length lower bound (default: 0)
  --length_upper_bound, -lu LENGTH_UPPER_BOUND
                        Length upper bound (default: 4294967296)
  --quality_threshold, -qt QUALITY_THRESHOLD
                        Length upper bound (default: 0)
  --logs_dir, -ld LOGS_DIR
                        Path to directory for logs (default: )
  --rewrite, -r         If should rewrite output (default: False)
```