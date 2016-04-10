#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>


void printMatrix(int p, int local_dim, double (matrix)[local_dim][local_dim][p][p]){
  int i;
  int j;
  int I;
  int J;

  for(I=0; I<p;I++){
        for(i=0; i<local_dim;i++){
            for(J=0;J<p;J++){
                for(j=0;j<local_dim;j++){
                    printf("%f \t",(matrix)[i][j][I][J]);
                }
            }
        printf("\n");
        }
    }
}

void printLocalMatrix(int local_dim, double matrix[local_dim][local_dim]){
    int i;
    int j;
    for(i=0;i<local_dim;i++){
        for(j=0;j<local_dim;j++){
            printf("%f \t",matrix[i][j]);
        }
        printf("\n");
    }
}

void printVector(int dim, double vector[dim]){
    int i;
    for(i=0;i<dim;i++){
        printf("%f \n",vector[i]);
    }
}

void matrixAdd(int dim, double matrixA[dim][dim], double matrixB[dim][dim], double result[dim][dim]){
    int i, j;
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            result[i][j]=matrixA[i][j]+matrixB[i][j];
        }
    }
}

void vectorSubtract(int dim, double vectorA[dim], double vectorB[dim], double result[dim]){
    int i;
    for(i=0;i<dim;i++){
        result[i]=vectorA[i]-vectorB[i];
    }
}

void matrixSubtract(int dim, double matrixA[dim][dim], double matrixB[dim][dim], double result[dim][dim]){
    int i, j;
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            result[i][j]=matrixA[i][j]-matrixB[i][j];
        }
    }
}

void vectorScalarMultiply(double scalar, int dim, double vector[dim], double result[dim]){
    int i;
    for(i=0;i<dim;i++){
        result[i]=scalar*vector[i];
    }
}

void matrixScalarMultiply(double scalar, int dim, double matrix[dim][dim], double result[dim][dim]){
    int i, j;
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            result[i][j]=scalar*matrix[i][j];
        }
    }
}

void matrixMultiply(int dim, double matrixA[dim][dim], double matrixB[dim][dim], double result[dim][dim]){
    int i, j, k;
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            for(k=0;k<dim;k++){
                result[i][j]+=(matrixA[i][k]*matrixB[k][j]);
            }
        }
    }
}

void randVector(int dim, double vector[dim]){
    srand((unsigned)time(NULL));
    int i;
    for(i=0; i<dim; i++){
        vector[i]=(double)rand();
    }
}

void randMatrix(int dim, double matrix[dim][dim]){
    srand((unsigned)time(NULL));
    int i,j;
    for(i=0; i<dim;i++){
        for(j=0; j<dim;j++){
            matrix[i][j]=(double)rand();
        }
    }
}

void onesVector(int dim, double vector[dim]){
    int i;
    for(i=0;i<dim;i++){
        vector[i]=1.0;
    }
}

void onesMatrix(int dim, double matrix[dim][dim]){
    int i,j;
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            matrix[i][j]=1.0;
        }
    }
}

double normVector(int dim, double vector[dim]){ // calculates euclidien vector norm
    int i;
    double norm;
    for(i=0;i<dim;i++){
        norm += vector[i]*vector[i];
    }
    norm = sqrt(norm);
    return norm;
}

double normMatrix(int dim, double matrix[dim][dim]){ // calculates L2,1 norm (sum of euclidien norms of columns)
    int i,j;
    double norm=0.0;
    double vectorNorm=0.0;
    for(j=0;j<dim;j++){
        vectorNorm=0.0;
        for(i=0;i<dim;i++){
            vectorNorm += (matrix[i][j]*matrix[i][j]);
        }
        vectorNorm = sqrt(vectorNorm);
        norm += vectorNorm;
    }
    return norm;
}

void copyMatrix(int dim, double matrixToCopy[dim][dim], double matrixToFill[dim][dim]){
    int i,j;
    for (i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            matrixToFill[i][j]=matrixToCopy[i][j];
        }
    }
}

void copyVector(int dim, double vectorToCopy[dim], double vectorToFill[dim]){
    int i;
    for(i=0;i<dim;i++){
        vectorToFill[i]=vectorToCopy[i];
    }
}
