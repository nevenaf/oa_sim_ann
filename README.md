# Simulated annealing for a construction of a random orthogonal array
### *** work in progress ***

Program in C for construction of random orthogonal arrays. It uses a simulated annealing technique and attempts to construct many arrays in parallel.

At the moment, the input is hard coded. Enter the number of columns (k) and the alphabet size (n) in run_oa_exchange.c. Also, initial temperature, cooling rate, and number of iterations can be adjusted here.

To compile run:
~~~
make
~~~

To execute:
~~~
./run_oa_exchange
~~~

There is also a possibility to extend an existing array. The **transpose** of the existing array should be entered in run_oa_exchange.c (elements of the first column, followed by the elements of the second column...).
