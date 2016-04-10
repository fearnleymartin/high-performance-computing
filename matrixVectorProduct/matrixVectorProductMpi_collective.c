// Matrix vector product in parallel
// Uses MPI collective functions such as scatter and gather
#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include "functions.c"

int main(int argc, char **argv){

  //MPI variables
  int rank; /* rank of the process */
  int size; /* number of processes */
  int source;
  int dest;
  int tag = 50;
  MPI_Status status;
  int root=0;

  //Matrix/vector variables
  int p = 5; //number of row bands
  int dim = 10;
  int local_dim=dim/p;

  //Initialise global matrix and vector and result
  double *matrix =NULL;
  double *vector=NULL;
  double *result=NULL;
  vector = malloc(dim*sizeof(double));
  result = malloc(dim*sizeof(double));

  //Initialise MPI
  MPI_Init (&argc, &argv); /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank); /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &size); /* get number of processes */

  if (rank==0){ //Process 0 initialises the matrix and vector

    //Initialise vector and matrix with random values
    matrix = malloc(dim*dim*sizeof(double));

    int i,j;
    srand((unsigned)time(NULL));
    for(i=0; i<dim;i++){
        vector[i]=(double)rand()/1000.0;
    }
    for (j=0; j<(dim*dim);j++){
        matrix[j] = (double)rand()/1000.0;
    }
    printf("Matrix: \n");
    printMatrix(dim, matrix); fflush(stdout);
    printf("Vector: \n");
    printVector(dim, vector);fflush(stdout);

  }

  //Initialise local_matrix and local_result
  double *local_matrix = malloc(dim*local_dim*sizeof(double));
  double *local_result = calloc(local_dim,local_dim*sizeof(double));


  MPI_Bcast(vector, dim, MPI_DOUBLE, root,MPI_COMM_WORLD);
  MPI_Scatter(matrix, (dim*local_dim), MPI_DOUBLE, local_matrix, (dim*local_dim), MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  //Local matrix multiplication
  int i,j;
  for(i=0;i<local_dim;i++){
    for(j=0;j<dim;j++){
      local_result[i] += (local_matrix[j+i*dim]*vector[j]);
    }
  }

  //Reassemble results from individual processes
  MPI_Gather(local_result, local_dim, MPI_DOUBLE, result, local_dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if(rank==0){
    printf("Result: \n");
    printVector(dim, result); fflush(stdout);
  }

  free(matrix);
  free(vector);
  free(local_matrix);
  free(local_result);

  MPI_Finalize();

}

