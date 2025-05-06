# ParaSeq: Parallel Sequence Alignment in Python
<small>
Università Politecnico di Milano<br>
Master's degree in Bioinformatics for Computational Genomics, Scientific Programming course
</small>

## Project
**Project proposal #6:** _Parallel implementation of global or local sequence alignment (Needlman-Wunsch or Smith-Waterman)_.<br>
This project provides a parallel implementation of the **Smith-Waterman** local sequence alignment algorithm in Python.

## How to run
**Clone repository** and move into the project folder:
```shell
git clone https://github.com/M1keCodingProjects/ParaSeq.git
cd ParaSeq
```

**Install dependencies** dependencies are inferred automatically from the pyproject.toml file, simply run this program:
```shell
pip install .
```
I recommend using **Python version 3.13.0** or higher.

**Run the program** and provide the necessary arguments:
```shell
python3 src\para_seq\main.py
```

## License
This project is licensed under the MIT License.

## Documentation
> Xia Z, et al. A Review of Parallel Implementations for the Smith-Waterman Algorithm. Interdiscip Sci. 2022;14(1):1-14. doi:10.1007/s12539-021-00473-0

_work in progress:_<br>
> The SW algorithm is used to find the best subsequence match between 2 sequences.