//Calculates page rank from initial matrix M
//Reads M from csv and exports result to csv
//Set N to number of pages on website (dimension of N) and p to the number of row bands (i.e. parallel computations)
//Currently p must divide N
//Example launch command: "mpiexec -n 10 pagerankmpi.exe" for example_test_matrix
//to launch with 10 processors
//Remember to specify the size of the entry matrix !
#include "functions1d.c"
#include <stdio.h>
#include <mpi.h>
#include <time.h>

char *input = "example_matrix.csv";
char *output = "pagerank.csv";

//Initialise with p processers
int main(int argc, char **argv){

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

    //Number of row bands
    int p = size;

    //Dimension of column of matrix M, i.e. number of web pages
    int N = 70;

    //Power iteration variables
    double *M; //adjacency matrix
    double d = 0.85; // damping factor
    double v_quadratic_error = 0.001; //quadratic error for v
    double *last_v; // last_v for starting loop and keeping track of previous v
    double *v; // vector used for power method
    double *E; // matrix of ones
    double *M_hat; //matrix used in power iteration method
    double *diff; // for storing v-last_v
    double *result; // for storing M_hat * v before rebuilding v

    //Local variables
    double *local_M_hat;
    int local_N;
    double *local_result;

    // Loop variables
    int loop = 1; //For syncronising loops between slave and master
    int numberOfIterations = 0; // For counting the number of iterations
    int i,j,k; //loop variables used in all processes

    //Initialise v because every process uses it
    v = malloc(N*sizeof(double)); //Allocate memory space for v

   // master process for initialising matrix and vector, and for running the loop and distributing the calculations
   if(rank==0){//Initialise problem in master process
        // Input of M
        M = malloc(N*N*sizeof(double)); //Allocate memory space for M
        importMatrix(input,N,M); //Affect coefficients to M from csv file
//        printf("%s : \n", "M");printMatrix(N,M);

        // Initialisation of a random v
        randVector(N,v); //Affect random coefficients to v
        vectorScalarMultiply(1.0/normVectorL1(N,v),N,v,v); // normalise v
//        printf("%s \n","Initial value of v"); printVector(N,v);

        //Initialisation of last_v
        last_v = malloc(N*sizeof(double)); //Allocate memory space for last_v
        onesVector(N, last_v); // Affect coefficients of last_v to 1
        vectorScalarMultiply(INFINITY, N, last_v, last_v); // Affect coefficients of last_v to infinity

        //Initialisation of E (matrix of ones)
        E = malloc(N*N*sizeof(double)); //Allocate memory space for E
        onesMatrix(N, E); //Affect coefficients of E to 1

        //Initialisation of M_hat
        M_hat = malloc(N*N*sizeof(double)); //Allocate memory space for M_hat
        //Calculation of M_hat
        matrixScalarMultiply(d,N,M,M); //Calculate M_hat step 1
        matrixScalarMultiply((1-d)/N,N,E,E); //Calculate M_hat step 2
        matrixAdd(N,M,E,M_hat); // Calculate M_hat

        // Once M_hat has been calculated, M and E are no longer needed in memory
        free(M); free(E);

        // Initialisation of result
        result = malloc(N*sizeof(double)); //Allocate memory space for result

        diff = malloc(N*sizeof(double)); //Allocate memory space to diff, used for calculating the difference between v and last_v

        printf("%s : \n", "M_hat");printMatrix(N,M_hat);fflush(stdout);

        // Send local_M to slaves (the matrices remain the same for every iteration)
        local_N=N/p; // get local_N
        for(dest=1; dest<size; dest++){ //Send to processes with rank != 0
            local_M_hat = M_hat + N*local_N*dest; // get pointer for local_M
            MPI_Send(&local_N,1,MPI_INT,dest,tag,MPI_COMM_WORLD); // send local_N to slaves
            MPI_Send(local_M_hat,N*local_N,MPI_DOUBLE,dest,2,MPI_COMM_WORLD); // send local_M to slaves
        }
        local_M_hat = malloc(N*local_N*sizeof(double));
        memcpy(local_M_hat,M_hat,N*local_N*sizeof(double)); //Get local M for processor 0
        free(M_hat); //Once local_M_hat s have been sent out, M_hat can be freed
        printf("%f \n",local_M_hat[0]);
    } //End if rank == 0

    else{ //rank != 0
        // Receives local_M from master (same for each iteration)
        MPI_Recv(&local_N,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
        local_M_hat = malloc(N*local_N*sizeof(double));
        MPI_Recv(local_M_hat,N*local_N,MPI_DOUBLE,0,2,MPI_COMM_WORLD,&status);
    } //end rank != 0

   //Start power iteration method
    while(loop==1){ //Loops until convergence is approached
        if(rank==0){ //Master for distributing tasks
            numberOfIterations++; //For seeing how many iterations are necessary for convergence
            printf("%s %d ##############################################\n","iteration", numberOfIterations);fflush(stdout);
            printf("%s \n","current value of v"); printVector(N,v);

            copyVector(N, v, last_v); // last_v  <- v
//            printf("%s \n","copied vector");fflush(stdout);
//
            //v = M_hat * v
//            Send local_v to slaves (different for each iteration) (local_M already sent)
            for(dest=1; dest<size; dest++){
                MPI_Send(v, N, MPI_DOUBLE,dest,3,MPI_COMM_WORLD); // send v to slaves
            }
//            printf("%s \n", "distributed v"); fflush(stdout);
//
            //Calculations and result for process 0
            for(i=0;i<local_N;i++){
                result[i]=0.0;
                for(j=0;j<N;j++){
                  result[i] += (local_M_hat[j+i*N] * v[j]); //Local matrix product
                }
            }
            // Receive local results and rebuild global result
            for(source=1; source<size; source++){//Local results for processes != 0
                MPI_Recv(&local_N, 1, MPI_INT, source, 4, MPI_COMM_WORLD, &status);
//                printf("%s \n", "local N received"); fflush(stdout);
                local_result = malloc(local_N*sizeof(double)); //Allocate memory space to local_result
                MPI_Recv(local_result, local_N, MPI_DOUBLE, source, 5, MPI_COMM_WORLD, &status);
//                printf("%s \n", "local result received"); fflush(stdout);
                for(j=0; j<local_N; j++){
                    result[j+local_N*(source)]=local_result[j]; //Building the global result from the local_results
                }
                free(local_result); //May be reconstructed with a different size
            }
//            printf("%s \n", "built global result"); fflush(stdout);
            copyVector(N, result, v); // v <- result
//            printf("%s \n","v has been updated with M_hat*v"); fflush(stdout);
            vectorSubtract(N,v,last_v,diff);

            if(numberOfIterations>150){ //For debugging purposes
                loop=0;
                MPI_Bcast(&loop,1,MPI_INT,0,MPI_COMM_WORLD);
            }

            //For communicating when to end the loop for the other
            if(normVectorL2(N,diff) < v_quadratic_error){//End loop when error is small enough
                loop=0;
                MPI_Bcast(&loop,1,MPI_INT,0,MPI_COMM_WORLD);
            }
            else{
                MPI_Bcast(&loop,1,MPI_INT,0,MPI_COMM_WORLD);
            }
        }//End rank ==0

        else{ //rank!=0 i.e.slaves
            // Receives v from master (different for each iteration)
            MPI_Recv(v,N,MPI_DOUBLE,0,3,MPI_COMM_WORLD,&status);
//            printf("%s \n", "received v");fflush(stdout);
            //Execute local calculation
            local_result = malloc(local_N*sizeof(double));

            for(i=0;i<local_N;i++){
                local_result[i]=0.0; //Initialise local_result with zeros
                for(j=0;j<N;j++){
                  local_result[i] += (local_M_hat[j+i*N] * v[j]); //Local matrix product
                }
            }
            MPI_Send(&local_N, 1, MPI_INT, 0, 4, MPI_COMM_WORLD);
//            printf("%s \n", "local N sent");fflush(stdout);
            MPI_Send(local_result, local_N, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);//Send local result to master
//            printf("%s \n", "local result sent");fflush(stdout);
            free(local_result);

            //For receiving instructions from process 0 on when to stop
            MPI_Bcast(&loop,1,MPI_INT,0,MPI_COMM_WORLD);
        } //End rank!=0
   } //End while loop

    if(rank==0){
        free(diff);
        free(last_v); // Can be freed here because only the master process uses it
        printf("%s \n", "v: "); printVector(N,v);
        exportVector(output,N,v); //Export result to csv file
    }
    else{
        free(local_M_hat); //After end of loop, local_m_hat is no longer needed
    }
    //Free up memory
    free(v);
    MPI_Finalize();
}
