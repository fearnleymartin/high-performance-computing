// matrix vector product parallel
#include <stdio.h>
#include <mpi.h>
#include <time.h>

void printMatrix(int *matrix, int dim){
  int i;
  int j;
  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
        printf("%d \t",matrix[i*dim+j]);
    }
    printf("\n");
  }
}

void printVector(int *vector, int dim){
    int i;
    for(i=0;i<dim;i++){
        printf("%d \n", vector[i]);
    }
}

int main(int argc, char *argv){

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
  int *matrix =NULL;
  int *vector=NULL;
  int *result=NULL;
  vector = malloc(dim*sizeof(int));
  result = malloc(dim*sizeof(int));

  //Initialise MPI
  MPI_Init (&argc, &argv); /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank); /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &size); /* get number of processes */

  if (rank==0){

    //Initialise vector and matrix with random values
    matrix = malloc(dim*dim*sizeof(int));

    int i;
    int j;

    srand(1);
    for(i=0; i<dim;i++){
        vector[i]=rand()/1000;
    }
    srand(2);
    for (j=0; j<(dim*dim);j++){
        matrix[j] = rand()/1000;
    }
    printf("Matrix: \n");
    printMatrix(matrix, dim);
    printf("Vector: \n");
    printVector(vector, dim);

  }

  //Initialise local_matrix and local_result
  int *local_matrix = malloc(dim*local_dim*sizeof(int));
  int *local_result = calloc(local_dim,local_dim*sizeof(int));


  MPI_Bcast(vector, dim, MPI_INT, root,MPI_COMM_WORLD);
  MPI_Scatter(matrix, (dim*local_dim), MPI_INT, local_matrix, (dim*local_dim), MPI_INT, root, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  int i;
  int j;
  for(i=0;i<local_dim;i++){
    for(j=0;j<dim;j++){
      local_result[i] += (local_matrix[j+i*dim]*vector[j]);
    }
  }

  MPI_Gather(local_result, local_dim, MPI_INT, result, local_dim, MPI_INT, 0, MPI_COMM_WORLD);
  if(rank==0){
    printf("Result: \n");
    printVector(result, dim);
  }

  free(matrix);
  free(vector);
  free(local_matrix);
  free(local_result);

  MPI_Finalize();

}

