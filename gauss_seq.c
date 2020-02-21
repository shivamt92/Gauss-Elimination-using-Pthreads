#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#define BILLION 1000000000L

void print(double **M,int r)
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

void printvect(double *f,int dim)
{   
    int i=0;
    for (i=0; i<dim; i++)
        printf(" %lf \n",f[i]);
}

void findswapmax(double **M, double *l,int c, int i,int r)
{       
    /* 
    M is the array of pointers to the A matrix
    l is the pointer to the array of b vector
    c is the column
    i is the row 
    r is the size of the matrix
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

    // Swapping 
    if (k!=i)
    {   
        double *m = M[k];

        M[k] = M[i];
        M[i] = m;

        //swapping b

        double temp = l[i];
        l[i] = l[k];
        l[k] = temp;

    }

}


int main(int argc, char *argv[])
{   
    struct timespec start,end;
    uint64_t diff;


    int r = atoi(argv[1]);
    double *A[r],mux,*oldA[r],*oldb;
    double *b,*res;               // CREATING AN ARRAY OF POINTERS
    int i,j,k,l;
    for (i=0; i<r; i++) 
    {
        A[i] = (double *)malloc(r * sizeof(double));           // ALLOCATING SPACE ON HEAP FOR THE MATRIX
        oldA[i] = (double *)malloc(r * sizeof(double));
    }

    b = (double *)malloc(r * sizeof(double));
    res = (double *)malloc(r * sizeof(double));
    oldb = (double *)malloc(r * sizeof(double));


     srand48(23);
     for (i=0; i<r; i++)           // Assiginig Values to A Matrix
    {
        for (j=0; j<r; j++)
        {   
            A[i][j] = drand48();                   
            oldA[i][j] = A[i][j];
        }
    }

    for (i=0; i<r; i++)             // ASSIGNING VALUES TO b vector
    {
    
        b[i] = drand48(); 
        oldb[i] = b[i];                  
        
    }

    // for (i=0; i<r; i++)
    // {
    //     for (j=0; j<r; j++)
    //     {   
    //         double x;
    //         scanf("%lf",&x);
    //         A[i][j] = x;                                  // ASSIGNING VALUES TO A Matrix
    //         oldA[i][j] = x;
    //     }
    // }

    // for (i=0; i<r; i++)
    // {
    //     double x;
    //     scanf("%lf",&x);
    //     b[i] = x;                   // ASSIGNING VALUES TO b vector
        
    // }

    clock_gettime(CLOCK_MONOTONIC,&start);
    // Starting Gaussian Elimination
    for (i=0; i<r-1; i++)
    {
        // printf("I was here");
        findswapmax(A,b,i,i,r);
        // printf("Through here \n");
        // printvect(b,r);
        // print(A,r);


        for(j=i+1; j<r; j++)
        {
            mux = A[j][i]/A[i][i];

            for(k=i; k<r; k++)
            {
                A[j][k] -=mux*A[i][k];
            }
            b[j] -= mux*b[i];
            // printvect(b,r);
        }

        // print(A,r);

    }

    // Starting Back Substitution


    for (i=r-1; i>-1; i-- )
    {
        
        res[i] = b[i]/A[i][i];


        for (j=0; j<i; j++)
        {
            b[j] -= A[j][i]*res[i]; 
        }
    }
    clock_gettime(CLOCK_MONOTONIC,&end);
    diff = BILLION * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;
    printf("\n TIME TAKEN FOR MATRIX MULTIPLICATION IS %lf  NANOSECONDS\n\n", (double) diff);
    
    double mul=0,error=0;
    int dim = r;
   for (i=0; i<r; i++)           // Calculating error
    {
        mul=0;
        for (j=0; j<dim; j++)
        {   
            mul+=oldA[i][j]*res[j];

        }
        error+= mul -oldb[i];

    }
    printf(" Error : %lf",error);

}



