# Experiments on Solving SVP on Rotations of Z^n

This repository contains the code and data for the experiments on the effects of the LLL and BKZ reduction algorithms on differently sampled bases of the integer matrix, as given in Section 6 of our paper "Just how hard are rotations of Z^n? Algorithms and cryptography with the simplest lattice". <insert link to paper>

The code used to implement the different basis generation algorithms (namely discrete gaussian-based sampling, unimodular matrix product sampling, and Bézout-coefficient-based sampling) are given in the file Z^n_sampling.py. The module makes use of NumPy and SageMath as dependencies, and should therefore be used in a SageMath virtual environment. 

The raw data for our experiements is given in the csv files. The naming convention for each column consists of some characters followed by a number: </br>
t - time taken (in seconds) </br>
max - largest matrix element </br>
lsq - shortest squared length of matrix </br>
1 - generate </br>
2 - LLL reduction </br>
k >= 3 - BKZ with block size k (applied successively)

For example, lsq1 is the shortest squared length of the generated matrix, max2 is the largest element of the LLL-reduced matrix and t5 is the time taken to run BKZ with block size 5 on the matrix.
