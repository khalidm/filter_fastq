# Overview

The filter_fastq tool removes reads containing a a given list of primer sequences.

# Licence

This program is released as open source software under the terms of [BSD-3-Clause License](https://raw.githubusercontent.com/bjpop/pdb_rename_chain/master/LICENSE).

# Installing

In the examples below, `$` indicates the command line prompt.

Install using pip and github.

```
pip install -U https://github.com/khalidm/filter_fastq/archive/master.zip

```

OR

Clone this repository:
```
$ git clone https://github.com/khalidm/filter_fastq
```

Move into the repository directory:
```
$ cd pdb_rename_chain
```

filter_fastq can be installed using `pip` in a variety of ways (`$` indicates the command line prompt):

1. Inside a virtual environment:
```
$ python3 -m venv filter_fastq
$ source filter_fastq/bin/activate
$ pip install -U /path/to/filter_fastq
```
2. Into the global package database for all users:
```
$ pip install -U /path/to/filter_fastq
```
3. Into the user package database (for the current user only):
```
$ pip install -U --user /path/to/filter_fastq
```

# General behaviour

Example usage:
```
filter_fastq -i <input_list_of_primers_sequences> -f1 <name_fastq_R1.fastq> -f2 <name_fastq_R2.fastq>

output:
name-1_fastq_R1.fastq
name-1_fastq_R2.fastq
```
