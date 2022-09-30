"""
Module containing implementations for sampling algorithms discussed in the paper.

This module makes use of SageMath. Thus, it has to be imported/run in a SageMath
virtual environment.
"""

import random
import numpy as np
import sage.all as sage

### HELPER FUNCTIONS ###

def xgcd2(a, b):
    """
    Runs the Extended Euclidean Algorithm on a and b, returning the gcd of the
    two numbers, d, followed by the Bézout coefficients x and y such that
    d = ax + by. Implementation adapted from that given by
    http://anh.cs.luc.edu/331/notes/xgcd.pdf

    Precondition: a and b are positive integers.
    """
    x_, x = 1, 0
    y_, y = 0, 1
    while b > 0:
        q = a // b
        x, x_ = x_ - q * x, x
        y, y_ = y_ - q * y, y
        a, b = b, a % b
    return a, x_, y_

def xgcd(L):
    """
    Runs the Extended Euclidean Algorithm on the integers in the list L,
    returning the gcd, d, of the integers of L = [i1, i2, ..., in], followed by
    the Bézout coefficients [x1, x2, ..., xn] such that
    d = x1 * i1 + x2 * i2 + ... + xn * in. Recursive with depth len(L).

    Precondition: L is a list of integers.
    """
    if len(L) == 1:
        return L[0], [1]
    d_, coefs = xgcd(L[:-1])
    d, x, y = xgcd2(d_, abs(L[-1]))
    return d, [x * coef for coef in coefs] + [np.sign(L[-1]) * y]

### SAMPLING ALGORITHMS ###

def discrete_gaussian(n, d, s, t):
    """
    Returns a sample of n d-dimensional row vectors (in the form of a
    nested list) whose entries are generated using a discrete gaussian
    distribution with parameter s. The algorithm uses rejection sampling within
    the integer range -t * s .. t * s.

    For the complete discrete gaussian based sampling discussed in the paper,
    the output of this algorithm (with d = n + 10) should be passed into the
    LLL-algorithm to get a square matrix representing a basis of the integer
    lattice.

    Precondition: n, d, s, t are positive integers.
    """
    A = [[0] * d for _ in range(n)]
    for i in range(n):
        for j in range(d):
            c = 0
            while not c:
                x = np.random.randint(-t * s, t * s + 1)
                c = np.random.binomial(1, np.e ** (- np.pi * x ** 2 / s ** 2))
            A[i][j] = x
    return A

def unimodular_product(n, B, L, d):
    """
    Samples and returns a n x n basis (in the form of a nested list) using
    the unimodular matrix product algorithm (basis sampling technique 2) given
    by the paper, where B is the size bound, L is the word length, and d is the
    size bound.

    Precondition: n, B, L and d are strictly positive integers and 2 <= d <= n.
    """
    M = sage.MatrixSpace(sage.IntegerRing(), n, sparse=True)
    A = M(sage.matrix.identity(n))
    for _ in range(L):
        cond = True
        while cond:
            U = np.random.randint(-B, B + 1, (d, d))
            cond = abs(np.linalg.det(U)) != 1
        y = np.zeros((n, n), dtype=object)
        K = random.sample(range(n), d)
        for i in range(d):
            for j in range(d):
                y[K[i], K[j]] = U[i, j] - (1 if i == j else 0)
        A += A * M(y)
    return np.array(A).tolist()

def bezout_sampling(n, B):
    """
    Samples and returns a n x n basis (in the form of a nested list) using the
    Bézout-coefficient-based construction (basis sampling technique 3) given by
    the paper, where B is the size bound.

    Precondition: n and B are strictly positive integers.
    """
    M = np.random.randint(-B, B + 1, (n, n))
    L = sage.matrix(M).adjugate().transpose()[0]
    d, coefs = xgcd(L)
    if d != 1:
        return bezout_sampling(n, B)
    A = np.vstack((M[1:], coefs)).tolist()
    return A
