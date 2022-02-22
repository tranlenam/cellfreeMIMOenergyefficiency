# On the total energy efficiency of cell-free massive MIMO

This repo contains the code for the algorithms presented in the following scientific paper:

Hien Quoc Ngo, Le-Nam Tran, Trung Q. Duong, Michail Matthaiou, Erik G. Larsson, "On the total energy efficiency of cell-free massive MIMO," IEEE Trans. Green Commun. and Network., vol. 2, no. 1, pp. 25-39, Mar. 2018.

## Instructions
The code makes use of [Yalmip](https://yalmip.github.io/) as a parser and [MOSEK](https://www.mosek.com/) as the internal convex conic solver for speed.
When running Algorithm 1 for many channel realizations and for large scale settings, you can replace the "inteferencevector" and "approxfunction" by their vectorised implementation, i.e., inteferencevectorvectorised" and "approxfunctionvectorised". You can refer to "Algorithm2.m" for how this can be done.

We also publish the code at [CodeOcean](https://codeocean.com/capsule/1803216/tree/v1) where you can run them online without the need to install the required packages.
