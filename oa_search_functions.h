#include<stdio.h>
#include<stdlib.h>

void print_int_matrix(int* M, int n, int r, FILE* fp);

void print_int_matrix_transpose(int* M, int n, int r, FILE* fp);

void make_random_object(int n, int k, int* object);
void make_random_population(int pop, int n, int k, int h, int* population, int* cost);

int sim_ann_v1(int n, int k, int* object, int old_k, int* old_object, int num_iterations, double initial_temp, double cool_rate);

void cost_vector(int n, int k, int* object, int* cv);

void num_missing(int n, int k, int* cv, int old_k, int* old_cv, int* missing);

void copy_array(int length, int* orig, int* copy);
