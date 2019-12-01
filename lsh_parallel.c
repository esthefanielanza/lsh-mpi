#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <mpi.h>
#include <stddef.h>

double RANDOM_ACCURACY = 0.7;
double THRESHOLD = 0.5;
long int LARGE_PRIME = 433494437;

typedef struct initialDataType {
  int nSets;
  int setSize;
  int stages;
  int buckets;
} initialDataType;

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

void generateRandomSets(int nSets, int setSize, int *sets) {
  srand(0);

  int i;

  for(i = 0; i < nSets * setSize; i++) {
    double sortedNumber = (double)rand() / (double)RAND_MAX;
    sets[i] = sortedNumber > RANDOM_ACCURACY;
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


void printVectorAsMatrix(int m, int n, int *vector) {
  int i;

  for(i = 0; i < m * n; i++) {
    printf("%d ", vector[i]);
    if(i != 0 && (i + 1) % n == 0) {
      printf("\n");
    }
  }
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


int h(int hashPosition, int setValue, int *coefs) {
  int start = 2 * hashPosition;
  return (int)((coefs[start] * (long) setValue + coefs[start + 1]) % LARGE_PRIME);
}

void calculateSignature(int *signature, int setSize, int signatureSize, int *set, int *coefs) {
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
void convertToSet(int* set, int *array, int start, int length) {
  int i, j = 0;

  for(i = start; i < start + length; i++) {
    if(array[i] == 1) {
      set[j] = i - start;
      j++;
    }
  }

  // Cleaning array
  while(j < length) {
    set[j] = 0;
    j++;
  }
}

void hashSignature(int *hashes, int currentHash, int *signature, int stages, int signatureSize, int buckets) {
  int rows = signatureSize / stages;

  for(int i = 0; i < signatureSize; i++) {
    int stage = min(i / rows, stages - 1);
    int start = stage + (currentHash * stages);
    hashes[start] = (int)((hashes[start] + (long) signature[i] * LARGE_PRIME) % buckets);
  }
}

void hashDataset(
  int *hashes,
  int setSize,
  int *sets,
  int stages,
  int buckets,
  int *coefs,
  int signatureSize,
  int partitionStart,
  int partitionEnd
) {
  int i;
  int *set = allocateVector(setSize);
  int *signature = allocateVector(signatureSize);

  for(i = partitionStart; i < partitionEnd; i++) {
    convertToSet(set, sets, i * setSize, setSize);
    calculateSignature(signature, setSize, signatureSize, set, coefs);
    hashSignature(hashes, i, signature, stages, signatureSize, buckets);
  }

  free(set);
  free(signature);
}

void printElementsPerBucket(int *hashes, int nSets, int stages, int buckets) {
  int i, j;
  int **counts = allocateMatrix(stages, buckets);

  for(i = 0; i < nSets; i++) {
    for(j = 0; j < stages; j++) {
      counts[j][hashes[i*stages + j]]++;
    }
  }

  printf("=== Buckets on each stage ===\n");
  printMatrix(stages, buckets, counts);

  deallocateMatrix(stages, counts);
}

MPI_Datatype createMpiInitialStruct() {
  int i;
  MPI_Datatype mpiInitialData;
  int blockLengths[4];
  MPI_Datatype types[4];
  MPI_Aint offsets[4];

  for(i = 0; i < 4; i++) {
    blockLengths[i] = 1;
    types[i] = MPI_INT;
  }

  offsets[0] = offsetof(initialDataType, nSets);
  offsets[1] = offsetof(initialDataType, setSize);
  offsets[2] = offsetof(initialDataType, stages);
  offsets[3] = offsetof(initialDataType, buckets);

  MPI_Type_create_struct(4, blockLengths, offsets, types, &mpiInitialData);
  MPI_Type_commit(&mpiInitialData);

  return mpiInitialData;
}

void setInitialDataOnProccess(initialDataType *initialData, double *start, MPI_Datatype mpiInitialData, int myRank, int nProcess) {
  int i;

  if(myRank == 0) {
    // Getting start time //
    *start = MPI_Wtime();

    // LSH Params //
    initialData->stages = 10;
    initialData->buckets = 10;

    // Dataset params ///
    initialData->nSets = 1000;
    initialData->setSize = 100;
    
    for(i = 1; i < nProcess; i++) {
      MPI_Send(initialData, 1, mpiInitialData, i, 0, MPI_COMM_WORLD);
    }
  } else {
    MPI_Recv(initialData, 1, mpiInitialData, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

int main () {
  // Parallel variables //
  int nProcess, myRank, i;
  double start, end;

  // Initialize the MPI enviroment //
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcess);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  initialDataType initialData;
  MPI_Datatype mpiInitialData = createMpiInitialStruct();
  setInitialDataOnProccess(&initialData, &start, mpiInitialData, myRank, nProcess);
  
  // Size of each partition //
  int partitionSize = nProcess > 0 ? ceil(initialData.nSets/nProcess) : initialData.nSets;

  int partitionSetSize = partitionSize * initialData.setSize;
  int *sets = allocateVector(partitionSetSize);
  if(myRank == 0) {
    int *allSets = allocateVector(initialData.nSets * initialData.setSize);

    // Generating entrace //
    generateRandomSets(initialData.nSets, initialData.setSize, allSets);

    // Master node should copy the hashes that it will work on //
    memcpy(sets, allSets, partitionSetSize * sizeof(int));

    // After that it should send to the other nodes the hashes that they need to work //
    for(i = 1; i < nProcess; i++) { 
      MPI_Send(&allSets[i * partitionSetSize], partitionSetSize, MPI_INT, i, 0, MPI_COMM_WORLD);
    }

    free(allSets);
  } else {
    // Slaves should receive their hashes and put that on the variable sets //
    MPI_Recv(sets, partitionSetSize, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  // Compute hash function coefficients //
  int signatureSize = getSignatureSize(initialData.stages);
  
  // Generating hash coefs //
  int *coefs = allocateVector(signatureSize*2);
  if(myRank == 0) {
    for(i = 0; i < signatureSize * 2; i++) {
      coefs[i] = (rand() % LARGE_PRIME) + 1;
    }
    
    // After generating coefs master node send all of them to the slaves //
    for(i = 1; i < nProcess; i++) {
      MPI_Send(coefs, signatureSize*2, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
  } else {
    // Slaves should reiceve coefs //
    MPI_Recv(coefs, signatureSize*2, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  

  // The size of the hash partition ( n of hashes per partition x hashSize (stages)) //
  int hashesSize = partitionSize * initialData.stages;
  
  // All nodes should handle the first N (partitionSize) hashes that they receive //
  int partitionStart = 0;
  int partitionEnd = partitionSize;

  if(myRank == 0) {
    // All hashes //
    int *hashes = allocateVector(initialData.nSets * initialData.stages);
    
    // Master calculate the hashes of their first N sets //
    hashDataset(
      hashes,
      initialData.setSize,
      sets,
      initialData.stages,
      initialData.buckets,
      coefs,
      signatureSize,
      partitionStart,
      partitionEnd
    );

    // Master waits until all proccess send their calculated hashes // 
    for(i = 1; i < nProcess; i++) {
      MPI_Recv(&hashes[hashesSize * i], hashesSize, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Print how the nodes where allocated between the buckets //
    printElementsPerBucket(hashes, initialData.nSets, initialData.stages, initialData.buckets);
  } else {
    // Just a part of the hashes //
    int *hashes = allocateVector(hashesSize);

    // Slaves calculate the hashes of their sets //
    hashDataset(
      hashes,
      initialData.setSize,
      sets,
      initialData.stages,
      initialData.buckets,
      coefs,
      signatureSize,
      0,
      partitionSize
    );

    // After it they send it to the master node //    
    MPI_Send(hashes, hashesSize, MPI_INT, 0, 0, MPI_COMM_WORLD);

    free(hashes);
  }

  if(myRank == 0) {
    end = MPI_Wtime();
    printf("\n%.6lf\n", (end - start)*1000.0);
  }

  MPI_Finalize();

  free(coefs);
  free(sets);
	return 0;
}