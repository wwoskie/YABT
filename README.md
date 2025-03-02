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