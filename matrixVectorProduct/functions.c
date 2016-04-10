#include <stdio.h>

void printMatrix(int dim, double *matrix){ //must be square matrix
  int i;
  int j;
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
        printf("%f \n", vector[i]);
    }
}
