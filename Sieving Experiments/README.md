# Sieving in $\mathbb{Z}^n$

## Downloading and Running Original Experiments

Clone this repository, navigate to the folder in which this README is in and run `make` to compile the executable. By running `./sieve-Zn [output_file]` you can run the same experiments we ran in our paper and store the outputs in [output_file]. Specifically, we run experiments to determine the effects the dimension $d$ of our integer lattice and the parameter of the discrete Gaussian $s$ we are sampling from have on the number of vectors the Gauss-Sieve algorithm samples, sieves, and compares and the time it takes the algorithm takes to run. We run experiments for every pair $(d,s)$, where $d \in [16,32,64]$ and $s \in [10,100,1000]$. 

## Running Other Experiments

If you would like to run other experiemnts — for different pairs $(d,s)$ — you can do so by changing the `s`  and `d` arrays in the `main` function.
