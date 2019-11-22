#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>

double RANDOM_ACCURACY = 0.7;
double THRESHOLD = 0.5;
long int LARGE_PRIME = 433494437;

//TODO: Criar desalocar matrix
int ** allocateMatrix(int nSets, int setSize) {
  int i;
  int **mat = (int **)calloc(nSets, sizeof(int*));

  for(i = 0; i < nSets; i++) {
    mat[i] = (int *)calloc(setSize, sizeof(int));
  }

  return mat;
}

void deallocateMatrix(int x, int **m) {
  int i;

  for(i = 0; i < x; i++) {
    free(m[i]);
  }

  free(m);
}

int *allocateVector(int size) {
  return (int *)calloc(size, sizeof(int));
}

void generateRandomSets(int nSets, int setSize, int **sets) {
  srand(0);

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
int getSignatureSize(int stages) {
  int r = (int) ceil(log(1.0 / stages) / log(THRESHOLD)) + 1;
  return r * stages;
}


int h(int hashPosition, int setValue, int **coefs) {
  return (int)((coefs[hashPosition][0] * (long) setValue + coefs[hashPosition][1]) % LARGE_PRIME);
}

void calculateSignature(int *signature, int setSize, int signatureSize, int *set, int ** coefs) {
  int i, j;

  for(i = 0; i < signatureSize; i++) {
    signature[i] = INT_MAX;
  }

  for(i = 0; i < setSize; i++) {
    for(j = 0; j < signatureSize; j++) {
      signature[j] = min(signature[j], h(j, set[i], coefs));
    }
  }
}

// Converts to an array where we count the 1's positions //
void convertToSet(int* set, int *array, int length) {
  int i, j = 0;

  for(i = 0; i < length; i++) {
    if(array[i] == 1) {
      set[j] = i;
      j++;
    }
  }

  // Cleaning array
  while(j < length) {
    set[j] = 0;
    j++;
  }
}

void hashSignature(int *hash, int *signature, int stages, int signatureSize, int buckets) {
  int rows = signatureSize / stages;

  for(int i = 0; i < signatureSize; i++) {
    int stage = min(i / rows, stages - 1);
    hash[stage] = (int)((hash[stage] + (long) signature[i] * LARGE_PRIME) % buckets);
  }
}

int ** hashDataset(int nSets, int setSize, int **sets, int stages, int buckets) {
  int i;
  int signatureSize = getSignatureSize(stages);

  printf("Signature size %d\n", signatureSize);

  // Compute hash function coefficients // 
  int **coefs = allocateMatrix(nSets, 2);
  int **hashes = allocateMatrix(nSets, stages);

  for(i = 0; i < signatureSize; i++) {
    coefs[i][0] = (rand() % LARGE_PRIME) + 1;
    coefs[i][1] = (rand() % LARGE_PRIME) + 1;
  }

  int *set = allocateVector(setSize);
  int *signature = allocateVector(signatureSize);

  printf("\n\n=== Computing coefficients ===\n");
  printMatrix(signatureSize, 2, coefs);

  printf("\n\n=== Calculating Hash ===\n");

  for(i = 0; i < nSets; i++) {
    convertToSet(set, sets[i], setSize);
    // printVector(setSize, set);

    calculateSignature(signature, setSize, signatureSize, set, coefs);
    // printVector(signatureSize, signature);

    hashSignature(hashes[i], signature, stages, signatureSize, buckets);
    printf("Hash[%d]:", i);
    printf(" : ");
    printVector(stages, hashes[i]);
  }

  free(set);
  free(signature);

  deallocateMatrix(nSets, coefs);
  return hashes;
}

void printElementsPerBucket(int **hashes, int nSets, int stages, int buckets) {
  int i, j;
  int **counts = allocateMatrix(stages, buckets);

  printf("\n\n=== Last stage position ===\n");
  for(i = 0; i < nSets; i++) {
    for(j = 0; j < stages; j++) {
      counts[j][hashes[i][j]]++;
      if(j == stages - 1) {
        printf("Set %d is on %d bucket\n", i, hashes[i][j]);
      }
    }
  }

  printf("\n\n=== Buckets on each stage ===\n");
  printMatrix(stages, buckets, counts);

  deallocateMatrix(stages, counts);
}

int main () {
  // Dataset params ///
  int nSets = 10;
  int setSize = 4;
  
  // Generating dataset // 
  int **sets = allocateMatrix(nSets, setSize);

  generateRandomSets(nSets, setSize, sets);

  printf("=== Sets ===\n");
  printMatrix(nSets, setSize, sets);
  
  // LSH Params //
  int stages = 2;
  int buckets = 4;

  // Generating hashes //
  int **hashes = hashDataset(nSets, setSize, sets, stages, buckets);

  printElementsPerBucket(hashes, nSets, stages, buckets);

  deallocateMatrix(nSets, sets);
  deallocateMatrix(nSets, hashes);
	return 0;
}