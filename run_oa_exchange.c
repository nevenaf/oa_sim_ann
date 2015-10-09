/*  *********COMPILATION*********

Method 1: make
          ./random_search
          
Method 2: gcc -std=c99 -fopenmp -c -o file_v1.o file_v1.c -lm  // for each suppor file
          gcc -std=c99 -fopenmp -o output_name file_main.c support_file1.o support_file2.o -lm 
          ./random_search
*/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<omp.h>
#include"oa_search_functions.h" 

#define num_iter_sa  5*1000000 //5*1000000 //4*1000000 //6*100000000; // 10^7
#define abs_initial_temp 0.05 //0.15 //original for PG(2,7) is 0.05
#define cool_rate 0.9999 //0.999 //0.99999 //original for PG(2,7) is 0.9999

#define num_tries 10
#define num_exchange 5


void print_combined_objects(int rows, int* mat1, int col1, int* mat2, int col2, FILE* f)
{
    int cols = col1+col2;
    int* new_mat = malloc(sizeof(int)*rows*cols);
    copy_array(rows*col1, mat1, new_mat);
    copy_array(rows*col2, mat2, new_mat+rows*col1);
    print_int_matrix(new_mat, rows, cols, f);
    free(new_mat);
}


int main(int argc, char **argv)
{

	int n=7;
	int k=8;
//	char* output_file_name=strcat(strcat("sim_ann_output_ca_k", k+'0'), strcat(strcat("_n", n+'0'), ".txt"));
	
	int old_k = 0;
	//char* input_file_name="input_array.txt";
	//int num_old_arrays = 1;
    
    
    // output file
	char output_file_name[256];
	sprintf(output_file_name, "sim_ann_output_ca_k%d_n%d.txt", k+old_k, n);
	
	printf("CA(k=%d, n=%d)\n", k+old_k,n);
	printf("-----------------\n\n");
    
	// parameters for sim. ann.
	int pop = pow(2,7); 
	
	printf("populaiton = %d, initial temp=%f \n", pop, abs_initial_temp);

	FILE* fp_res = fopen(output_file_name, "a+");
    
    int* population; // population of matrices of size h x r
    int* cost; // cost of each member of the population
    int min_cost;
    int* best_sol;
    int* best_sol_cost=malloc(sizeof(int)*pop);
    double initial_temp;
 
 	int h = n*n;
 	int total = k*(k-1)/2 + old_k*k;
    
    population = malloc(sizeof(int)*h*(k)*pop);
    cost = malloc(sizeof(int)*pop); 
    best_sol = malloc(sizeof(int)*h*(k)*pop);
      
	int save_idx;
    
    srand(time('\0')); // seed timer with the current time
    min_cost = h*total;
    
    int* old_object;
    if (old_k > 0)
    {
        old_object = malloc(sizeof(int)*h*old_k);    
        for(int i=0; i<n; i++)
        {
            for(int j=0; j<n; j++)
            {
                old_object[i*n+j] = i;
                old_object[h+i*n+j] = j;
            }
        }
        fprintf(fp_res, "\nexisting (old) object: \n");
        print_int_matrix(old_object, h, old_k, fp_res);
    }


    for(int try=0; try<num_tries && min_cost>0 ; try++)
    {
        // look for the new best solution
        for(int p=0; p<pop; p++) best_sol_cost[p]=h*total;
        initial_temp = abs_initial_temp;

        make_random_population(pop, n,k,h, population,cost);
        min_cost = h*total;
        
        for(int exch = 0 ; exch<num_exchange && min_cost>0; exch++)
        {

            if(min_cost>0)
            {
                // run sim. ann. on each member of population
                int p;
                //if(exch>0) initial_temp=initial_temp*3/4;
                #pragma omp parallel private(p) shared(n,k,population, old_k, old_object,  initial_temp, cost, min_cost)
                {
                    #pragma omp for
                    for (p=0; p<pop; p++)
                    {
                        cost[p]= sim_ann_v1( n, k, population+h*k*p, old_k, old_object,  num_iter_sa,initial_temp, cool_rate);
                        if(min_cost>cost[p]){ min_cost = cost[p];}
                    }
                } 
                printf("finished exchange exch=%d min_cost=%d \n", exch, min_cost);
                printf("some costs: ");
                for(p=0; p<15; p++){printf("%d ", cost[p]);}
                printf("\n");
                /*
                for(int p=0; p<pop; p++)
                {
                    if(best_sol_cost[p]>cost[p]) 
                    {
                        copy(h*r, population+h*r*p, best_sol+h*r*p);
                        best_sol_cost[p]=cost[p];
                        if(min_cost>  best_sol_cost[p]) min_cost = best_sol_cost[p];  
                    }
                }
                fprintf(fp, "\nexchange=%d Sim. ann. done! cost of the first %d: ", exch, pop);
                for(int p=0; p<pop; p++) fprintf(fp, "%d, ", cost[p]);
                fprintf(fp, "\n");
                */

            }
        }

        //if(min_cost==0)
        {
            fprintf(fp_res, "n=%d k=%d \n", n,k);
            
            save_idx=0;
            for(int p=1; p<pop; p++)
            {
                if(cost[save_idx]>cost[p]) save_idx=p;
            } 
               
            fprintf(fp_res, "\nobject: \n");
            
            print_combined_objects(h, old_object, old_k, population+h*(k)*save_idx, k, fp_res);
            //print_int_matrix(population+h*(k)*save_idx, h,(k), fp_res);
            
            fprintf(fp_res, "\n cost= %d \n", cost[save_idx]);
            fprintf(fp_res, "\n");
            fflush(fp_res);
            
        }
	}

	printf("Results saved in file: %s \n", output_file_name);
	//fclose(fp);
	fclose(fp_res);
	
	free(population);
	free(cost);

	free(best_sol);
	free(best_sol_cost);
	
	if(old_k > 0)
	{
	    free(old_object);
   }

	return min_cost;

}

