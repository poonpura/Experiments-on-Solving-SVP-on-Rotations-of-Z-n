# Basis Reduction Experiments

This repository contains the code and data for the experiments on the effects of the LLL and BKZ reduction algorithms on differently sampled bases of the integer matrix, as given in Section 6 of our paper "Just how hard are rotations of Z^n? Algorithms and cryptography with the simplest lattice". <insert link to paper>

The code used to implement the different basis generation algorithms (namely discrete Gaussian-based sampling, unimodular matrix product sampling, and Bézout-coefficient-based sampling) are given in the file Z^n_sampling.py. The module makes use of NumPy and SageMath as dependencies, and should therefore be used in a SageMath virtual environment. 

The raw data for our experiements is given in the csv files. The naming convention for each column consists of some characters followed by a number: </br>
t - time taken (in seconds) </br>
max - largest absolute value of a matrix element </br>
lsq - shortest squared length of a vector in the matrix </br>
1 - the initial generation algorithm </br>
2 - LLL reduction </br>
k >= 3 - BKZ with block size k (applied successively)

For example, lsq1 is the shortest squared length of a vector in the generated matrix, max2 is the largest absolute value of an element of the LLL-reduced matrix, and t5 is the time taken to run BKZ with block size 5 on the matrix.

## How to use sampling algorithms

The sampling algorithms used in our experiments are provided in the `Zn_sampling.py` file, with documentation for the specific sampling functions provided as docstrings. Be sure to run the module in a SageMath virtual environment. To make use of them, import the file as a module (the import may take a couple of minutes).
```
>>> import Zn_sampling as zns
>>> zns.discrete_gaussian(3, 3, 10, 10)
[[-4, -2, 0], [-6, -4, 2], [-1, -2, 4]]
```
To view documentation in shell, run 
```
>>> help(zns)
```
