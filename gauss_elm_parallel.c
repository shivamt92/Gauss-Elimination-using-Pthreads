/*  ------------------- SHIVAM TEWARI - NETID: st949 -------------------------------------------
    ------------------- ECE 5720: INTRODUCTION TO PARALLEL COMPUTING ---------------------------
    ----------------- ASSIGNMENT 2 : GAUSS ELIMINATION USING PTHREADS --------------------------
    ---------------- COMPILE: gcc gauss_elm_parallel.c -pthreads -o gpar -----------------------
    -------------------------- ARG 1: SIZE OF THE ARRAY ----------------------------------------
    -------------------------- ARG 2: NUMBER OF PTHREADS ---------------------------------------
*/

#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#define BILLION 1000000000L


struct gauss_t                      // Global Variable thread of struct gauss_t
{
    int col;
    int row;
    double **m;
    double *b;
    double *res;
    int dim;
    int numthread;
    
}thread;

void* Gauss_elm(void* arg)              // Gauss - Elimination
{      
    int* actual;
    actual = (int *) arg;
    double **M = thread.m;
    double *b = thread.b;
    int row = thread.row;
    int dim = thread.dim;
    int numthreads = thread.numthread;
    int tid = *actual;

    int i,j,k;
    double mux;
    
    for (k = row + 1+tid; k<dim; k+=numthreads)
    {
    mux = M[k][row]/M[row][row];
    
        for(i=row; i<dim; i++)
        {
        M[k][i] -=mux*M[row][i];
        }
    b[k] -= mux*b[row];
    }
}


void* Back_sub(void* arg)               // Back Substitution
{
    int *actual;
    actual = (int *) arg;
    double **M = thread.m;
    double *b = thread.b;
    double *res = thread.res;
    int row = thread.row;
    int dim = thread.dim;
    int numth = thread.numthread;
    int tid = *actual;
    int i;
    for (i = 0+tid ; i<row; i+=numth)
    {
        b[i]-= M[i][row]*res[row];
    }
   
}


void print(double **M,int r)            // Function to print Matrix
{
    int i,j;
    for (i=0; i<r; i++)
    {
        for (j=0; j<r; j++)
        {   
            printf("%lf  ",M[i][j]);                 // ASSIGNING RANDOM VALUES TO BOTH THE MATRICES
        }
        printf("\n");
    }
    return ;
}


void findswapmax(double **M, double *l,int c, int i,int r)       // Function for partial pivoting
{       
    /* 
    M : array of pointers to the A matrix
    l : pointer to the array of b vector
    c : column
    i : row 
    r : size of the matrix
    */
    int j,k;
    double max = fabs(M[i][c]);
    k=i;

    for(j =i+1; j<r; j++)
    {   
        if (fabs(M[j][c]) > max)
        {
            k = j;
            max = fabs(M[j][c]);
        }
    }
    if (k!=i)                  // Swapping 
    {   
        double *m = M[k];

        M[k] = M[i];
        M[i] = m;

        double temp = l[i];    // Swapping b
        l[i] = l[k];
        l[k] = temp;
    }
}

void initializematrix(double **m, double **oldm,int d)              // Initialize a matrix
{   
    int i;
    for (i=0; i<d; i++)           // ALLOCATING SPACE ON HEAP FOR A MATRIX
    {
        m[i] = (double *)malloc(d * sizeof(double));
        oldm[i] = (double *)malloc(d * sizeof(double));           
    }
}

void assignvaluesmatrix(double **m,double **oldm,int d)             // Assign Values manually to MAtrix
{   
    int i,j;  
    for (i=0; i<d; i++)           // Assiginig Values to A Matrix
    {
        for (j=0; j<d; j++)
        {   
            double x;
            scanf("%lf",&x);
            m[i][j] = x;                   
            oldm[i][j] = x;
        }
    }
}

void random_assign_matrix(double **m,double **oldm,int d)           // Assign value randomly to matrix
{   
    int i,j;  
    srand48(22);
     for (i=0; i<d; i++)           // Assiginig Values to A Matrix
    {
        for (j=0; j<d; j++)
        {   
            m[i][j] = drand48();                   
            oldm[i][j] = m[i][j];
        }
    }
    
}

void random_assign_b(double *a, double *olda,int d)                 // Assign values randomly to b vector
{   
    int i;
    srand48(23);
    for (i=0; i<d; i++)             // ASSIGNING VALUES TO b vector
    {
        a[i] = drand48();                   
        olda[i] = a[i];
    }
}

void initialize_vectors(double *a, double *olda, double *re, int d)     
{
    a = (double *) malloc(d * sizeof(double));
    re = (double *) malloc(d * sizeof(double));
    olda = (double *) malloc(d * sizeof(double));

}

void assign_to_b(double* b, double* oldb,int dim)       // ASSIGNING VALUES TO b vector manually
{
    int i;
    for (i=0; i<dim; i++)             
    {
        double x;
        scanf("%lf",&x);
        b[i] = x;                   
    }

}
int main(int argc, char **argv)
{   
    struct timespec start,end;
    uint64_t diff;

    int num_th = atoi(argv[2]);                     
    int dim = atoi(argv[1]);
    double *A[dim],*oldA[dim], *b,*oldb, *res;
    int i,j,k,l,m,n;

    for (i=0; i<dim; i++)                               // ALLOCATING SPACE ON HEAP FOR A MATRIX
    {
        A[i] = (double *)malloc(dim * sizeof(double));
        oldA[i] = (double *)malloc(dim * sizeof(double));           
    }

    b = (double *) malloc(dim * sizeof(double));            // Allocating space for b vector and resultant vector
    oldb = (double *) malloc(dim * sizeof(double)); 
    res = (double *) malloc(dim * sizeof(double));

    // initializematrix(A,oldA,dim);
    // initialize_vectors(b,oldb,res,dim);

    // assignvaluesmatrix(A,oldA,dim);
    // assign_to_vec(b,oldb,dim);
    
    random_assign_matrix(A,oldA,dim);
    random_assign_b(b,oldb,dim);

    pthread_t thread_id[num_th];                    // Initializing Threads
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

    int r =dim;
    thread.m = A;
    thread.b = b;
    thread.dim = dim;
    thread.numthread = num_th;
    thread.res = res;

    int index[num_th];
    for (i=0; i<num_th; i++)
    {
        index[i] = i;
    }
    printf("\n ------------------------------------ ECE 5720: ASSIGNMENT 2 -----------------------------------------------------------\n");
    printf("\n ------------------------------------ FINDING SOLUTIONS TO Ax=b USING GAUSS ELIMINATION AND BACK SUBSTITUTION ----------\n ");
    printf("\n ------------------------------------ REDUCING A MATRIX INTO ROW ECHELON FORM USING PTHREADS ----------------------------\n");
    clock_gettime(CLOCK_MONOTONIC,&start);

    for (i=0; i<dim-1; i++)                       // Gauss Elimination using pthreads
    {
        findswapmax(A,b,i,i,r);                  // Finding the max and pivoting
        
        thread.row = i;

            for (k=0; k<num_th; k++)
            {   
            pthread_create(&thread_id[k],&attr,Gauss_elm,(void *) &index[k]);       
            }
            for (m=0; m<num_th; m++)
            {
                int rc = pthread_join(thread_id[m],NULL);
                if (rc){
                    printf("ERORO %d %d",rc,m);
                    exit(-1);
                }
            }
    }
    
    printf("\n ------------------------------------ FINISHED GAUSSIAN ELIMINATION ----------------------------------------------------\n");

    printf("\n------------------------------------- BACKSUBSTITUTING USING PTHREADS --------------------------------------------------\n");

    
    for (i=dim-1; i>-1; i-- )           // Backsubstitution using Pthreads
    {
        res[i] = b[i]/A[i][i];
        thread.row = i;

        for (k=0; k<num_th; k++)
        {   
            pthread_create(&thread_id[k],&attr,Back_sub,(void *) &index[k]);       
        }
        for (m=0; m<num_th; m++)
        {
        int rc = pthread_join(thread_id[m],NULL);
        }
    }
    
            
    clock_gettime(CLOCK_MONOTONIC,&end);
    diff = BILLION * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;
    printf("\n TIME TAKEN FOR FINDING SOLUTIONS TO X IS %lf  NANOSECONDS\n\n", (double) diff);

    // Calculate Error in parallel

    double error = 0,mul=0;

    for (i=0; i<dim; i++)                 // Calculating error
    {
        mul=0;
        for (j=0; j<dim; j++)
        {   
            mul+=oldA[i][j]*res[j];
        }
        error+= mul -oldb[i];

    }
    printf(" Error : %lf \n\n",error);

}