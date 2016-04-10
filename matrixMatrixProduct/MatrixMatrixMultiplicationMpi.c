// matrix matrix product in parallel
// master divides up matrix into p**2 sub matrixes and distributes each block multiplication to a processor
// each local block matrix is only stored in its corresponding processor. When another block matrix is needed, the mpi message system is used to send it this block matrix when needed
// launch command example "mpiexec -n 4 matrixmatrixmultiplicationmpi.exe"
// where 4 should be replaced by the number of processors (i.e. the number of sub matrices you divide your matrix into)
// Note: p must divide dim, i.e. the square root of the number of processors should divide the global matrix dimension
#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include "functions.c"

//Matrix/vector variables
int p; //number of row bands
int dim; //Matrix dimension
int local_dim; //Dimension of block matrices


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
  int p = 2; //number of row bands and should be equal to the sqrt of the number of processes
  int dim = 4;
  int local_dim=dim/p;

  //Initialise global matrix and vector and result
  double matrixA[local_dim][local_dim][p][p];
  double matrixB[local_dim][local_dim][p][p];
  double result[local_dim][local_dim][p][p];

  double localMatrixA[local_dim][local_dim];
  double localMatrixB[local_dim][local_dim];
  double localResult[local_dim][local_dim];

  double localMatrixATemp[local_dim][local_dim];
  double localMatrixBTemp[local_dim][local_dim];

  //Initialise localResult
  int i,j;
  for(i=0;i<local_dim;i++){
    for(j=0;j<local_dim;j++){
        localResult[i][j]=0.0;
    }
  }

  //Initialise MPI
  MPI_Init (&argc, &argv); /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank); /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &size); /* get number of processes */

  // Create communicators
  MPI_Comm row_comm;
  MPI_Comm column_comm;
  int I, J;
  I=(rank)%p;
  J=(rank)/p;
  MPI_Comm_split(MPI_COMM_WORLD,I,J,&row_comm);
  MPI_Comm_split(MPI_COMM_WORLD,J,I,&column_comm);
//  int row_rank, row_size;
//MPI_Comm_rank(row_comm, &row_rank);
//MPI_Comm_size(row_comm, &row_size);
//
//printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n", rank, size, row_rank, row_size);




  if (rank==0){

    //Initialise vector and matrix with random values
    int i,j,I,J;

    srand(time(NULL));
    for(I=0; I<p;I++){
        for(i=0; i<local_dim;i++){
            for(J=0;J<p;J++){
                for(j=0;j<local_dim;j++){
                    matrixB[i][j][I][J]=(double)rand()/1000.0;
                }
            }
        }
    }
    for(I=0; I<p;I++){
        for(i=0; i<local_dim;i++){
            for(J=0;J<p;J++){
                for(j=0;j<local_dim;j++){
                    matrixA[i][j][I][J]=(double)rand()/1000.0;
                }
            }
        }
    }

    printf("Matrix A: \n");
    printMatrix(p, local_dim, matrixA);
    printf("Matrix B: \n");
    printMatrix(p, local_dim, matrixB);

    // Send submatrixes AI,J BI,J CI,J to processors
    for(dest=1;dest<size;dest++){
        // The processes are numbered starting from 1 (0 is the master process)
        I=(dest)%p;
        J=(dest)/p;
        for(i=0;i<local_dim;i++){
            for(j=0;j<local_dim;j++){
                localMatrixA[i][j]=matrixA[i][j][I][J];
                localMatrixB[i][j]=matrixB[i][j][I][J];
            }
        }
        MPI_Send(localMatrixA, local_dim*local_dim, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
        MPI_Send(localMatrixB, local_dim*local_dim, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
    }

    // Get A1,1
    for(i=0;i<local_dim;i++){
        for(j=0;j<local_dim;j++){
            localMatrixA[i][j]=matrixA[i][j][0][0];
            localMatrixB[i][j]=matrixB[i][j][0][0];
        }
    }
  }
  else{
    // Receive submatrixes AI,J BI,J CI,J
    source=0;
    MPI_Recv(localMatrixA, local_dim*local_dim, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
    MPI_Recv(localMatrixB, local_dim*local_dim, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
  }

    // Do block calculations
    int K;
    for(K=0;K<p;K++){
        if(K==J){
            memcpy(localMatrixATemp, localMatrixA, local_dim*local_dim*sizeof(double));
        }
        MPI_Bcast(localMatrixATemp,local_dim*local_dim, MPI_DOUBLE, K, row_comm);
        if(K==I){
            memcpy(localMatrixBTemp, localMatrixB, local_dim*local_dim*sizeof(double));
        }
        MPI_Bcast(localMatrixBTemp,local_dim*local_dim, MPI_DOUBLE, K, column_comm);
        //A = A + Btemp*Ctemp
        matrixMultiply(local_dim, localMatrixATemp, localMatrixBTemp, localResult);
    }

    // Send local results back to master
    if(rank!=0){
        dest=0;
        MPI_Send(localResult, local_dim*local_dim, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
    }
    else{ // Rank==0
        //Get A0,0
        int i,j;
        for(i=0;i<local_dim;i++){
            for(j=0;j<local_dim;j++){
                result[i][j][0][0]=localResult[i][j];
            }
        }
        //Get local_results
        for(source=1;source<size;source++){
            MPI_Recv(localResult, local_dim*local_dim, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
            I=(source)%p;
            J=(source)/p;
            int i,j;
            for(i=0;i<local_dim;i++){
                for(j=0;j<local_dim;j++){
                    result[i][j][I][J]=localResult[i][j];
                }
            }
        }
        printf("Result: \n");
        printMatrix(p, local_dim, result);
    }

  // MPI finalize
  MPI_Comm_free(&row_comm);
  MPI_Comm_free(&column_comm);
  MPI_Finalize();
}

