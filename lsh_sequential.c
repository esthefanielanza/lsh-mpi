#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

double RANDOM_ACCURACY = 0.7;
double THRESHOLD = 0.5;
long int LARGE_PRIME = 433494437;

//TODO: Criar desalocar matrix
int ** allocateMatrix(int nSets, int setSize) {
  int i;
  int **mat = (int **)malloc(nSets * sizeof(int*));

  for(i = 0; i < nSets; i++) {
    mat[i] = (int *)malloc(setSize * sizeof(int));
  }

  return mat;
}

int *allocateVector(int size) {
  return (int *)malloc(size *sizeof(int));
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

void printVector(int size, int *vector) {
  int i;

  for(i = 0; i < size; i++) {
    printf("%d ", vector[i]);
  }

  printf("\n");
}

int min(int x, int y) { 
  return (x < y) ? x : y;
}

int * hash(int stages, int buckets, int *set, int setSize) {
  int i;
  int *hash = allocateVector(stages);
  int rows = setSize / stages;

  for(i = 0; i < setSize; i++) {
    int stage = min(i / rows, stages - 1);
    hash[stage] = (int)(( hash[stage] + (long) set[i] * LARGE_PRIME ) % buckets);
  }

  return hash;
}

// Compute the size of signature //
int getSignatureSize(stages, setSize) {
  int r = (int) ceil(log(1.0 / stages) / log(THRESHOLD)) + 1;
  return r * stages;
}

void hashDataset(int **hashes, int nSets, int setSize, int **sets, int stages, int buckets) {
  int i;
  int signatureSize = getSignatureSize(stages, setSize);

  printf("Signature size %d\n", signatureSize);

  // Compute hash function coefficients // 
  int **coefs = allocateMatrix(nSets, 2);

  for(i = 0; i < setSize; i++) {
    coefs[i][0] = rand() + 1;
    coefs[i][1] = rand() + 1;
  }

  printMatrix(nSets, 2, coefs);
}



int main () {
  // Dataset params ///
  int nSets = 10;
  int setSize = 10;
  
  // Generating dataset // 
  int **sets = allocateMatrix(nSets, setSize);

  generateRandomSets(nSets, setSize, sets);

  // printMatrix(nSets, setSize, sets);
  
  // LSH Params //
  int stages = 2;
  int buckets = 5;

  // Generating hashes //
  int *hashes = allocateVector(nSets);
  hashDataset(&hashes, nSets, setSize, sets, stages, buckets);

	return 0;
}