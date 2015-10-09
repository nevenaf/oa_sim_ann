/**
    My original implementation. 
    It does take into account that it may be given some arrays.
    Annealing: done with the worst column choice.
**/

#include"oa_search_functions.h"
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<omp.h>

void print_int_matrix(int* M, int n, int r, FILE* fp)
{
    for(int row=0; row < n; row++)
    {
        fprintf(fp, "[");
        for(int col=0; col<r-1; col++)
        {
            fprintf(fp, "%d, ", *(M+row+col*n));
        }
        fprintf(fp, "%d], \n", *(M+(r-1)*n+row));
    } 
}

void print_int_matrix_transpose(int* M, int n, int r, FILE* fp)
{
    int h = n*n;
    for(int col=0; col<r; col++)
    {
        fprintf(fp, "[");
        for(int row=0; row < n-1; row++)
        {
            fprintf(fp, "%d, ", *(M+row+col*n));
        }
        fprintf(fp, "%d], \n", *(M+col*h+(n-1)));
    } 
}

// Returns a random integer i <= uniform(i,m) <= m 
unsigned uniform(unsigned i,unsigned m)
{
    //return (unsigned)(floor(((double)rand()/(1.0+RAND_MAX))*(m+1-i)+i));
    return (unsigned)(rand()%(m+1-i)+i);
} 

// Cost function: binary vector of all pairs covered by each pair of columns
void cost_vector(int n, int k, int* object, int* cv)
{
	int total = (k)*(k-1)/2;
	int count = 0;
	int h = n*n;
	int a, b;
	for(int i=0; i<h*total; i++) cv[i]=0;
	
	for(int c1=0; c1<k; c1++)
	{
		for(int c2=c1+1; c2<k; c2++)
		{
			for(int l=0; l<h; l++) 
			{
				a = object[c1*h+l];
				b = object[c2*h+l];
				cv[count*h+a*n+b]+=1;
			}
			count ++;	
		}
	}
}

// Cost function: binary vector of all pairs covered by a pair of columns
// such that the first colum is from the old array "old_object" and the second
// colum is from the second array called "object"
void old_cost_vector(int n, int k, int* object, int old_k, int* old_object, int* old_cv)
{
    int count = 0;
	int h = n*n;
	int a, b;
	for(int i=0; i<h*k*old_k; i++) old_cv[i]=0;
	
	for(int c1=0; c1<old_k; c1++)
	{
		for(int c2=0; c2<k; c2++)
		{
			for(int l=0; l<h; l++) 
			{
				a = old_object[c1*h+l];
				b = object[c2*h+l];
				old_cv[count*h+a*n+b]+=1;
			}
			count ++;
		}
	}
}

void num_missing(int n, int k, int* cv, int old_k, int* old_cv, int* missing)
{
	int covered = 0; 
	int count = 0;
	int h = n*n;
	for(int c=0; c<k+old_k; c++) missing[c]=0;
	for(int c1 = 0 ; c1<k; c1++)
	{
		for(int c2=c1+1; c2<k; c2++)
		{
			covered = 0;
			for(int l=0; l<h; l++)
			{
				if(cv[count*h+l]>0)
				{
					covered += 1;
				}
			}
			missing[old_k+c1]+=n*n-covered;
			missing[old_k+c2]+=n*n-covered;
			count ++;
		}
		for(int c2=0; c2<old_k; c2++)
		{
			covered = 0;
			for(int l=0; l<h; l++) if(old_cv[(c2*k+c1)*h+l]>0) covered ++;
			missing[old_k+c1]+=n*n-covered;
			missing[c2]+=n*n-covered;
		}
	}
}

// make a random object of size n^2 x (k), so that each colum is a random LS
// n = number of symbols
// k = number of squares
// object = address of final object *** don't forget to deallocate it ***
// object is an array of length (k) * n^2 (entered by columns)
void make_random_object(int n, int k, int* object)
{    
    int h = n*n;
    // enter elements in all other columns
    for(int c=0; c<k; c++)
    {
    	for(int s=0; s<n; s++)
		{
			for(int l=0; l<n; l++)
			{
				object[c*h+s*n+l] = s;
			}
		}
    }
    
    // randomly rearrange elements in each side of the object
    int i;
    int j;
    int swap;
    for(int s=0; s<k; s++)
    {
        for (i=0; i<h; i++) 
        {
            j = uniform(i,h-1);
            swap = object[s*h+i];
            object[s*h+i] = object[s*h+j];
            object[s*h+j] = swap;
        }
    }

}

void make_random_population(int pop, int n, int k, int h, int* population, int* cost)
{
    int p;
    #pragma omp parallel private(p) shared(pop,n,k,h,population,cost)
    {

        #pragma omp for
        for (p=0; p<pop; p++)
        {
            make_random_object(n, k, population+p*h*k);
            cost[p]=0;
        }
    }
}

int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

// ascending order
void sort(int h, int* cv, int* cv_enum)
{
    for(int l=0; l<h; l++)
    {
        cv_enum[2*l]=cv[l];
        cv_enum[2*l+1]=l;
    }
    qsort(cv_enum, h, 2*sizeof(int), cmpfunc);
}

void copy_array(int length, int* orig, int* copy)
{
	for(int i=0; i<length; i++) copy[i]=orig[i];
}


// pick col with priority reverse proportional to cost
int pick_col(int len_forbid, int len_in_use, int* cost, int cv_min, int cv_max)
{
    int rn=-1;
    while(rn==-1)
    {
        rn = rand()%len_in_use;
        if((double)rand()/ ((double)RAND_MAX+1) > (1.0-(cv_max-cost[len_forbid+rn])/(cv_max-cv_min+0.5))) rn=-1;
    }
    return rn;
}

int cost_remove(int row, int col, int n, int k, int* obj, int* cv, int* idx_mat)
{
	int h = n*n;
	int val, idx;
	int cost = 0;
	for(int c=0; c<k; c++)
	{
		if(c!=col)
		{
			if(c<col) val = obj[c*h+row]*n+obj[col*h+row];
			else val = obj[col*h+row]*n+obj[c*h+row];
		}
		idx = idx_mat[c*(k)+col];
		if(cv[idx*h+val]==1) cost += 1;
	} 

	return cost;
}

int swap_cost(int row1, int row2, int col, int n, int k, int* obj, int* cv, int* idx_mat, int old_k, int* old_obj, int* old_cv)
{
	int h = n*n;
	int val11, val12, val21, val22, idx;
	int cost = 0;
	
	// check if we are swapping different values
	if(obj[col*h+row1]!=obj[col*h+row2])
	{
		for(int c=0; c<k; c++)
		{
			// check if the values in c are the same in row1 and row2
			if((c!=col) && (obj[c*h+row1]!=obj[c*h+row2]))
			{
			
				if(c<col) 
				{ 
					val11 = obj[c*h+row1]*n+obj[col*h+row1]; 
					val12 = obj[c*h+row2]*n+obj[col*h+row2];
					val21 = obj[c*h+row1]*n+obj[col*h+row2];
					val22 = obj[c*h+row2]*n+obj[col*h+row1];
				}	
				else 
				{
					val11 = obj[col*h+row1]*n+obj[c*h+row1]; 
					val12 = obj[col*h+row2]*n+obj[c*h+row2];
					val21 = obj[col*h+row1]*n+obj[c*h+row2];
					val22 = obj[col*h+row2]*n+obj[c*h+row1];
				}
			
				idx = idx_mat[c*(k)+col];
				if(cv[idx*h+val11]==1) cost += 1;
				if(cv[idx*h+val12]==1) cost += 1;
				if(cv[idx*h+val21]==0) cost -= 1;
				if(cv[idx*h+val22]==0) cost -= 1;
			}
		}
		
		// repeat the same for the old_object
		for(int c=0; c<old_k; c++)
		{
			// check if the values in c are the same in row1 and row2
			if(old_obj[c*h+row1]!=old_obj[c*h+row2])
			{
				val11 = old_obj[c*h+row1]*n+obj[col*h+row1]; 
				val12 = old_obj[c*h+row2]*n+obj[col*h+row2];
				val21 = old_obj[c*h+row1]*n+obj[col*h+row2];
				val22 = old_obj[c*h+row2]*n+obj[col*h+row1];
			
				idx = c*k+col;
				if(old_cv[idx*h+val11]==1) cost += 1;
				if(old_cv[idx*h+val12]==1) cost += 1;
				if(old_cv[idx*h+val21]==0) cost -= 1;
				if(old_cv[idx*h+val22]==0) cost -= 1;
			}
		}
		 
	}
	return cost;
}

void old_remove_value(int col, int row, int n, int k, int* obj, int old_k, int* old_obj, int* old_cv, int* missing)
{
	int h = n*n;
	int val, idx;

	for(int c=0; c<old_k; c++)
	{
		val = old_obj[c*h+row]*n+obj[col*h+row];

		idx = c*k+col;
		if(old_cv[idx*h+val]==1) 
		{
			missing[old_k+col] += 1;
			missing[c] += 1;
		}
		old_cv[idx*h+val] -= 1;
	} 
}


void remove_value(int col, int row, int n, int k, int* obj, int old_k, int* cv, int* missing, int* idx_mat)
{
	int h = n*n;
	int val, idx;

	for(int c=0; c<k; c++)
	{
		if(c!=col)
		{
			if(c<col) val = obj[c*h+row]*n+obj[col*h+row];
			else val = obj[col*h+row]*n+obj[c*h+row];
		
			idx = idx_mat[c*(k)+col];
			if(cv[idx*h+val]==1) 
			{
				missing[old_k + c] += 1;
				missing[old_k + col] += 1;
			}
			cv[idx*h+val] -= 1;
		}
	} 
}

void add_value(int new_val, int col, int row, int n, int k, int* obj, int old_k, int* cv, int* missing, int* idx_mat)
{
	int h = n*n;
	int val, idx;
	for(int c=0; c<k; c++)
	{
		if(c!=col)
		{
			if(c<col) val = obj[c*h+row]*n+new_val;
			else val = new_val*n+obj[c*h+row];
		
			idx = idx_mat[c*(k)+col];
		
			if(cv[idx*h+val]==0) 
			{
				missing[old_k + c] -= 1;
				missing[old_k + col] -= 1;
			}
			cv[idx*h+val] += 1;
		}
	} 
}

void old_add_value(int new_val, int col, int row, int n, int k, int* obj, int old_k, int* old_obj, int* old_cv, int* missing)
{
	int h = n*n;
	int val, idx;
	for(int c=0; c<old_k; c++)
	{
		val = old_obj[c*h+row]*n+new_val;
		idx = c*k+col;

		if(old_cv[idx*h+val]==0) 
		{
			missing[old_k + col] -= 1;
			missing[c] -= 1;
		}
		old_cv[idx*h+val] += 1;
	} 
}

void make_adjustments(int col, int row1, int row2, int n, int k, int* object, int* cv, int old_k, int* old_object, int* old_cv, int* missing, int* idx_mat)
{
    int h = n*n;
	remove_value(col, row1, n, k, object, old_k, cv, missing, idx_mat);
	remove_value(col, row2, n, k, object, old_k, cv, missing, idx_mat);
	old_remove_value(col, row1, n, k, object, old_k, old_object, old_cv, missing);
	old_remove_value(col, row2, n, k, object, old_k, old_object, old_cv, missing);
	add_value(object[col*h+row2], col, row1, n, k, object, old_k, cv, missing, idx_mat);
	add_value(object[col*h+row1], col, row2, n, k, object, old_k, cv, missing, idx_mat);
	old_add_value(object[col*h+row2], col, row1, n, k, object, old_k, old_object, old_cv, missing);
	old_add_value(object[col*h+row1], col, row2, n, k, object, old_k, old_object, old_cv, missing);
	
}

double acceptance_prob(int cost,  double temp)
{
    return exp(-(double)cost)*temp;
}

int sim_ann_v1(int n, int k, int* object, int old_k, int* old_object, int num_iterations, double initial_temp, double cool_rate)
{
	int h = n*n;
	// copy of the best found object
	int* best_object = malloc(sizeof(int)*h*(k));
	int min_found_cost=0;
	
  	// evaluate the existing object
    int total = (k)*(k-1)/2;
    int* cv = malloc(sizeof(int)*total*h);
    int* missing = malloc(sizeof(int)*(k+old_k));
    
    // addition: cost vector with old_object
    int* old_cv;
    if (old_k > 0)
    {
        old_cv = malloc(sizeof(int)*old_k*k*h);
        // addition: initialize the old_cv, and adjust num_missing
        old_cost_vector(n, k, object, old_k, old_object, old_cv);
    }
    
    cost_vector(n, k, object, cv);
    num_missing(n, k, cv, old_k, old_cv, missing);
    
	copy_array(h*(k), object, best_object);
	for(int c=0; c<k+old_k; c++) min_found_cost += missing[c];
	min_found_cost = min_found_cost/2;
//	printf("missing: ");
//    print_int_matrix(missing, 1, k, stdout);
//	printf("min found cost=%d \n", min_found_cost);
//	printf("old cv: \n");
//	print_int_matrix(old_cv, h, k*old_k, stdout);
//	printf("cv: \n");
//	print_int_matrix(cv, h, total, stdout);
	
	int min_cost;
	int max_cost;
	int col, row11, row12, row21, row22, row1, row2;
	int total_cost1, total_cost2, total_cost;
	int swap;
	int curr_cost = min_found_cost;
	double current_temp = initial_temp;
	int good;

	// build index matrix 
	int* idx_mat = malloc(sizeof(int)*(k)*(k));
	int count = 0;
	for(int c1=0; c1<k; c1++)
	{
		idx_mat[c1*(k)+c1] = -1;
		for(int c2=c1+1; c2<k; c2++)
		{	
			idx_mat[c1*(k)+c2] = count;
			idx_mat[c2*(k)+c1] = count;
			count ++;
		}
	}
	long unsigned iter=0;
    while(iter<num_iterations && min_found_cost>0) //(curr_cost > 0) 
    {
        // evaluate the current object and find min and max costs
        min_cost = h*(k+old_k);
        max_cost = 0;
        for(int c=0; c<k; c++)
        {
            if(min_cost>missing[old_k+c]) min_cost=missing[old_k+c];
            if(max_cost<missing[old_k+c]) max_cost=missing[old_k+c];
        }
        
        
        // pick a random (with priority) column
        col = pick_col(old_k, k, missing, min_cost, max_cost);

		row11 = rand()%h;
		row12 = row11;
		while(row12 == row11) row12 = rand()%h;
		
		row21 = rand()%h;
		row22 = row21;
		while(row22 == row21) row22 = rand()%h;
		
		//printf("column = %d \n", col);
		//printf("first choice of rows: %d, %d \n", row11, row12);	
		//printf("second choice of rows: %d, %d \n", row21, row22);	
		total_cost1 = swap_cost(row11, row12, col, n, k, object, cv, idx_mat, old_k, old_object, old_cv);
		total_cost2 = swap_cost(row21, row22, col, n, k, object, cv, idx_mat, old_k, old_object, old_cv);
		//printf("total_cost1 = %d, total_cost2=%d \n", total_cost1, total_cost2);
        

        if(total_cost1<total_cost2)
        {
            total_cost = total_cost1;
			row1 = row11;
			row2 = row12;
        }
        else
        {
            total_cost = total_cost2;
            row1 = row21;
			row2 = row22;
        }
        
        
        good=0;
        if(total_cost <= 0) good=1;
        else if(acceptance_prob(total_cost,  current_temp)> ( ((double)rand())/ ((double)RAND_MAX+1) )) 
        {
            good=1; 
            current_temp = current_temp*cool_rate;
            if(current_temp < 0.0001) current_temp=0.0001;
        }
        
        
        if(good==1)
        {
        
        	// make changes to cv and missing
            make_adjustments(col, row1, row2, n,k, object, cv, old_k, old_object, old_cv, missing, idx_mat);
            
            // accept swap
            swap = object[col*h+row1];
            object[row1+col*h] = object[row2+col*h];
            object[row2+col*h] = swap;
            curr_cost += total_cost;

            if(min_found_cost>curr_cost)
            {
                min_found_cost=curr_cost;
                copy_array(h*(k), object, best_object);
                //printf("new_cost=%d  \n", curr_cost);
            }
            //printf("total_cost = %d new cost = %d iter=%lu \n", total_cost, curr_cost, iter);
            //printf("current array: \n");
            //print_int_matrix(object, h, k, stdout);
            //printf("cv: \n");
            //print_int_matrix(cv, h, total, stdout);
            //printf("missing: ");
            //print_int_matrix(missing, 1, k+old_k, stdout);
            //printf("\n\n");
        }
        
        iter++;
        //if(iter%1000==0) printf("iter=%lu \n", iter);
    }
	


    free(idx_mat);
    free(cv);
    if (old_k > 0)
    {
        free(old_cv);
    }
    free(missing);
    copy_array(h*(k), best_object, object);
    free(best_object);
    
    return min_found_cost;
}


