// matrix vector product parallel
#include <stdio.h>
#include <mpi.h>
#include <time.h>

void printMatrix(int *matrix, int dim){ //must be square matrix
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
//Initialise MPI
  MPI_Init (&argc, &argv); /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank); /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &size); /* get number of processes */
  //Matrix/vector variables
  int p = size-1; //number of row bands, each processor treats one band and the final processor is the master
  int dim = 10; //dimension of matrix row/vector
  //Initialise global matrix and vector and result
  int *matrix =NULL;
  int *vector=NULL;
  int *result=NULL;
  vector = malloc(dim*sizeof(int));
  //Initialise local_matrix and local_result
  int local_dim;
  int *local_result=NULL;
  int *local_matrix=NULL;

  if (rank==0){//The master processor. Splits matrix, distributes calculations, and reassembles matrix at the end

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
    //Initialise result
    result = malloc(dim*sizeof(int));

    // Create and send vector, local_dim and local matrices to each processor
    for(dest=1; dest<size; dest++){
        //Calculate local dim
        if (dim % p ==0){ // the number of row band is a divisor of the dimension of the matrix
            local_dim=dim/p;
            //Get local matrices
            local_matrix = matrix + dim*local_dim*(dest-1);
        }
        else{ // the number of row bands is not a divisor of the dimension of the matrix
            if(dest != (size-1)){ //All the row bands except the last
              local_dim=dim/p;
              //Get local matrices
              local_matrix = matrix + dim*local_dim*(dest-1);
            }
            else{ // The last row band is the size needed to complete the matrix
                local_dim=(dim-((dim/p)*(p-1)));
                //Get local matrix
                local_matrix = matrix + dim*(dim/p)*(dest-1);
            }
        }

        //Send data
        MPI_Send(&local_dim,1,MPI_INT,dest,tag,MPI_COMM_WORLD);
        MPI_Send(vector,dim,MPI_INT,dest,tag,MPI_COMM_WORLD);
        MPI_Send(local_matrix,dim*local_dim,MPI_INT,dest,tag,MPI_COMM_WORLD);

    }

    //Receive the local_dim and local_results from processors
    for(source=1; source<size; source++){
        MPI_Recv(&local_dim, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
        int *local_result = malloc(local_dim*sizeof(int));
        MPI_Recv(local_result, local_dim, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
        // Construction of result from local results
        if(source != size-1){ // All the local results apart from the last one, which may be of a different size
            for(j=0; j<local_dim; j++){
                result[j+local_dim*(source-1)]=local_result[j];
            }
        }
        else if(source==size-1){ // For the band, the size may be different, so a different indexing is used
                //((dim/p)*(p-1)) takes us to the index where the final band starts and we fill up from there
            for(j=0; j<local_dim; j++){
                result[j+((dim/p)*(p-1))]=local_result[j];
            }
        }
    }
    // Print final result
    printf("Result: \n");
    printVector(result, dim);
  }

  else { //i.e. rank != 0. Each processor calculate a local part of the product
    //Receive vector, local_dim and local matrix
    MPI_Recv(&local_dim, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
    int *local_matrix = malloc(dim*local_dim*sizeof(int));
    MPI_Recv(vector, dim, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
    MPI_Recv(local_matrix, dim*local_dim, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
    //Initialise local result with 0s
    int *local_result = malloc(local_dim*sizeof(int));
    int k;
    for(k=0;k<local_dim;k++){
      local_result[k]=0;
    }
    //Each processor does its local calculation
    int i;
    int j;
    for(i=0;i<local_dim;i++){
        for(j=0;j<dim;j++){
          local_result[i] += (local_matrix[j+i*dim] * vector[j]);
        }
    }
    //Each processor sends its local_dim and its local result to the master processor
    MPI_Send(&local_dim, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
    MPI_Send(local_result, local_dim, MPI_INT, 0, tag, MPI_COMM_WORLD);
  }

  // Free variables and finalise MPI
  free(matrix);
  free(vector);
  free(local_matrix);
  free(local_result);
  free(result);
  MPI_Finalize();
}

