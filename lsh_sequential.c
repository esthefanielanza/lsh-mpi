#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

double RANDOM_ACCURACY = 0.7;

//TODO: Criar desalocar matrix
int ** allocateMatrix(int nSets, int setSize) {
  int i;
  int **mat = (int **)malloc(nSets * sizeof(int*));

  for(i = 0; i < nSets; i++) {
    mat[i] = (int *)malloc(setSize * sizeof(int));
  }

  return mat;
}

void generateRandomSets(int nSets, int setSize, int **sets) {
  srand(time(NULL));

  int i, j;

  for(i = 0; i < nSets; i++) {
    for(j = 0; j < setSize; j++) {
      double sortedNumber = (double)rand() / (double)RAND_MAX;
      sets[i][j] = sortedNumber > RANDOM_ACCURACY;
    }
  }
}

void printMatrix(int nSets, int setSize, int **sets) {
  int i, j;

  for(i = 0; i < nSets; i++) {
    for(j = 0; j < setSize; j++) {
      printf("%d ", sets[i][j]);
    }
    printf("\n");
  }
}

int main () {
  int nSets = 10;
  int setSize = 10;
  
  int **sets = allocateMatrix(nSets, setSize);

  generateRandomSets(nSets, setSize, sets);

  printMatrix(nSets, setSize, sets);
  
	return 0;
}