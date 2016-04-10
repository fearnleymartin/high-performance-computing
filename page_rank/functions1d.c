#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "csvparser.h"
#include "csvparser.c"
#include "csvwriter.h"
#include "csvwriter.c"

void printMatrix(int dim, double *matrix){ //must be square matrix
  int i,j;
  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
        printf("%f ",matrix[i*dim+j]);
    }
    printf("\n");
  }
}

void printVector(int dim, double *vector){
    int i;
    for(i=0;i<dim;i++){
        printf("%f \n",vector[i]);
    }
}

void matrixAdd(int dim, double *matrixA, double *matrixB, double *result){
    int i,j;
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            result[i*dim+j]=matrixA[i*dim+j]+matrixB[i*dim+j];
        }
    }
}

void vectorSubtract(int dim, double *vectorA, double *vectorB, double *result){
    int i;
    for(i=0;i<dim;i++){
        result[i]=vectorA[i]-vectorB[i];
    }
}

void matrixSubtract(int dim, double *matrixA, double *matrixB, double *result){
    int i, j;
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            result[i*dim+j]=matrixA[i*dim+j]-matrixB[i*dim+j];
        }
    }
}

void vectorScalarMultiply(double scalar, int dim, double *vector, double *result){
    int i;
    for(i=0;i<dim;i++){
        result[i]=scalar*vector[i];
    }
}

void matrixScalarMultiply(double scalar, int dim, double *matrix, double *result){
    int i,j;
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            result[i*dim+j]=scalar*matrix[i*dim+j];
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

void randVector(int dim, double *vector){
    srand((unsigned)time(NULL));
    int i;
    for(i=0; i<dim; i++){
        vector[i]=((double)rand())/((double)RAND_MAX);
    }
}

void randMatrix(int dim, double *matrix){
    srand((unsigned)time(NULL));
    int i;
    for(i=0; i<dim*dim;i++){
        matrix[i]=((double)rand())/((double)RAND_MAX);
    }
}

void onesVector(int dim, double *vector){
    int i;
    for(i=0;i<dim;i++){
        vector[i]=1.0;
    }
}

void onesMatrix(int dim, double *matrix){
    int i;
    for(i=0;i<dim*dim;i++){
        matrix[i]=1.0;
    }
}

double normVectorL2(int dim, double *vector){ // calculates euclidien vector norm (L2)
    int i;
    double norm=0.0;
    for(i=0;i<dim;i++){
        norm += vector[i]*vector[i];
    }
    norm = sqrt(norm);
    return norm;
}

double normVectorL1(int dim, double *vector){ // calculates L1 vector norm
    int i;
    double norm=0.0;
    for(i=0;i<dim;i++){
        norm += vector[i];
    }
    return norm;
}


double normMatrix(int dim, double *matrix){ // calculates L2,1 norm (sum of euclidien norms of columns)
    int i,j;
    double norm=0.0;
    double vectorNorm=0.0;
    for(j=0;j<dim;j++){
        vectorNorm=0.0;
        for(i=0;i<dim;i++){
            vectorNorm += (matrix[i*dim+j]*matrix[i*dim+j]);
        }
        vectorNorm = sqrt(vectorNorm);
        norm += vectorNorm;
    }
    return norm;
}

void copyMatrix(int dim, double *matrixToCopy, double *matrixToFill){
    int i,j;
    for (i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            matrixToFill[i*dim+j]=matrixToCopy[i*dim+j];
        }
    }
}

void copyVector(int dim, double *vectorToCopy, double *vectorToFill){
    int i;
    for(i=0;i<dim;i++){
        vectorToFill[i]=vectorToCopy[i];
    }
}

int importMatrix(char *input, int dim, double *m){
    int i = 0;
    int j =  0;
    //                                   file, delimiter, first_line_is_header?
    CsvParser *csvparser = CsvParser_new(input, ";", 0);
    CsvRow *row;

    while ((row = CsvParser_getRow(csvparser)) ) {
//    	printf("==NEW LINE==\n");
        const char **rowFields = CsvParser_getFields(row);
        for (j = 0 ; j < CsvParser_getNumFields(row) ; j++) {
//            printf("FIELD: %s\n", rowFields[j]);
            m[i*dim+j]=atof(rowFields[j]);
        }
//		printf("\n");
        CsvParser_destroy_row(row);
        i++;
    }
    CsvParser_destroy(csvparser);

    return 0;

}

int exportVector(char *output, int dim, double *vector){

	CsvWriter *csvWriter = CsvWriter_new(output, ",", 0);
	int i;
	for (i = 0 ; i < dim ; i++) {
        char s[50];
        sprintf(s,"%f",vector[i]);
        if (CsvWriter_writeField(csvWriter, s)) {
            printf("Error: %s\n", CsvWriter_getErrorMessage(csvWriter));
            return 1;
        }
		CsvWriter_nextRow(csvWriter);
	}
	CsvWriter_destroy(csvWriter);

	return 0;
}

// int main(int argc, char **argv){
////     double *M;
////     double *v;
////     int N  = 70;
////     v = malloc(N*sizeof(double));
////     randVector(N,v);
////     exportVector("pageRank.csv",N,v);
//
////     M = malloc(N*N*sizeof(double));
////     importMatrix("finance_utile_M.csv",N,M);
////     printMatrix(N,M);
////    int matrixA[2][2] = {{1,2},{3,4}};
////    int matrixB[2][2] = {{5,10},{15,20}};
////    double matrixA[2][2];
////    randMatrix(2,matrixA);
////    double matrixB[2][2];
////    onesMatrix(2,matrixB);
////
////    double vector[2];
////    onesVector(2,vector);
////    printf("%f ",normVector(2,vector));
////
////    double scalar=2;
//////    printLocalMatrix(2,matrixA);
////    printLocalMatrix(2,matrixB);
////    double newMatrix[2][2];
//
////    matrixSubtract(2, matrixA, matrixB, newMatrix);
////    matrixScalarMultiply(scalar, 2, matrixA, newMatrix);
////    printLocalMatrix(2,newMatrix);
////    printf("%f",normMatrix(2,matrixB));
//
//}
