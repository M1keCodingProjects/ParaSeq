# ParaSeq: Parallel Sequence Alignment in Python
<small>
Università Politecnico di Milano<br>
Master's degree in Bioinformatics for Computational Genomics, Scientific Programming course
</small>

## Project
**Project proposal #6:** _Parallel implementation of global or local sequence alignment (Needlman-Wunsch or Smith-Waterman)_.<br>
This project provides a parallel implementation of the **Smith-Waterman** local sequence alignment algorithm in Python.
**Author:** _Michele Ferrari_

## How to run
**Clone repository** and move into the project folder:
```shell
git clone https://github.com/M1keCodingProjects/ParaSeq.git
cd ParaSeq
```

Make sure you are using **Python version 3.13.0** or higher, to check run:
```shell
python --version
```

**Install dependencies**, which are inferred automatically from the pyproject.toml file.
Simply run this command:
```shell
pip install .
```

**Run the program** as a module and provide the necessary arguments:
```shell
python -m src.para_seq.main ACGTGTG .\data\good.fasta -m 2 -mm 1 -g 2 -qp 1
```

**Read the full output** a summary output will be shown in the terminal, but if you want
to read the full list of alignments head on over to the ```.\output``` folder. Read
the documentation to see how to **change the ouput path**.
If you executed the example above you should find results that match what you'll find in
```.\output\example_output.txt```, perhaps in a different order.

## Documentation
The tool will accept **direct DNA sequences** or **FASTA** (.fasta, .fa) **file paths**,
but in both cases only A, C, G, T and N are allowed as nucleotides. When using FASTA file
paths it's possible to specify additional arguments (```-tp``` and ```-qp```) respectively
allowing the user to pick the target and query sequence from a specific position in the
respective files.
If a single file path is used both positions will be used to index from the same file, but
keep in mind that the tool **doesn't allow the alignment of identical sequences** (even when
passed as raw DNA strings).

The match, mismatch and gap (constant, applied equally to gap creation and extension) scores
are passed to the tool as **non-negative integers**, then the tool will interpret the mismatch
and gap scores as penalties, and use them to subtract from the alignment score.

The optional ```-o``` argument allows the user to specify a file path for the output file
of the tool, containing all the local alignments that were found and the alignment score.
If the argument is not specified the file will be available in the ```.\output\``` subfolder.

To further control the output the user can employ the ```-ma``` and ```-ls``` optional
arguments to limit the amount of alignments shown in the terminal summary output and to
limit the lengths of the shown aligned sequences.

This tool achieves the **parallelization** of the Smith-Waterman algorithm in two distinct
steps of the pipeline:
- Matrix filling step: as detailed in:
    > Xia Z, et al. A Review of Parallel Implementations for the Smith-Waterman Algorithm. Interdiscip Sci. 2022;14(1):1-14. doi:10.1007/s12539-021-00473-0

    the matrix filling step can be parallelized by iterating through the matrix not by rows or
    columns but by **antidiagonals**. Since each cell needs the cell directly above, the one
    to its left and the one diagonally above and to the left to compute its score, each cell
    in an antidiagonal can be evaluated independently. As such, each cell in an antidiagonal
    is sent to a different MultiProcessing Process (the specifics of dispatch are handled
    automatically with a **MultiProcessing Pool**, using as many cores as there are available). 

- Backtracking step: when reconstructing the optimal local alignments the starting point is
    the cell with the highest score. In case of ties optimal alignments can start from all
    of the tied cells with a maximum score and the reconstruction of one path is never
    affected by the others, which means that each group of alignments starting from the same
    position can be computed **in parallel**. Once again, each group is dispatched to a
    MultiProcessing Process via a **MultiProcessing Pool**.

**Interesting detail #1:** in order for the directions matrix to occupy less space in memory
I introduced an optimization where all the directions in a cell are represented by different
bits in a 3-bit binary flag. Essentially, each cell in the matrix is a **uint8** where the
3 least significant bits are set as such:
- b001 : can go up from here
- b010 : can go diagonally up-left from here
- b100 : can go left from here

Therefore for example a cell containing both the "up" and "diagonal" directions would have
the value b00000011 (3, in decimal). This greatly reduces the memory footprint of the
directions matrix from a more naive implementation with a list of numbers or strings.

**Interesting detail #2:** the alignments for 2 runs on identical data will be the same but
most likely in a different order. This is because the alignments are collected in a set in
order to **eliminate true duplicates** (exactly the same aligned subsequences, from exactly the
same starting positions).

## License
This project is licensed under the MIT License.