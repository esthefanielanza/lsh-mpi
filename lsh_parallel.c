#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <mpi.h>

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


int h(int hashPosition, int setValue, int *coefs, int signatureSize) {
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
      signature[j] = min(signature[j], h(j, set[i], coefs, signatureSize));
    }
  }
}

// Converts to an array where we count the 1's positions //
void convertToSet(int* set, int *array, int start, int length) {
  int i, j = 0;

  for(i = start; i < start + length; i++) {
    // printf("%d ", array[i]);
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

  // printf(" = ");
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
  int nSets,
  int setSize,
  int *sets,
  int stages,
  int buckets,
  int *coefs,
  int signatureSize,
  int partitionStart,
  int partitionEnd
) {
  int i, j;
  int *set = allocateVector(setSize);
  int *signature = allocateVector(signatureSize);

  // printf("\n\n=== Computing coefficients ===\n");
  // printMatrix(signatureSize, 2, coefs);

  // printf("\n\n=== Calculating Hash ===\n");

  for(i = partitionStart; i < partitionEnd; i++) {
    convertToSet(set, sets, i * setSize, setSize);
    // printVector(setSize, set);
    calculateSignature(signature, setSize, signatureSize, set, coefs);
    // printVector(signatureSize, signature);
    hashSignature(hashes, i, signature, stages, signatureSize, buckets);
  }

  // printVectorAsMatrix(partitionEnd - partitionStart, stages, hashes);
  free(set);
  free(signature);
}

void printElementsPerBucket(int *hashes, int nSets, int stages, int buckets) {
  int i, j;
  int **counts = allocateMatrix(stages, buckets);

  // printf("\n\n=== Last stage position ===\n");
  for(i = 0; i < nSets; i++) {
    for(j = 0; j < stages; j++) {
      counts[j][hashes[i*stages + j]]++;
      // if(j == stages - 1) {
      //   printf("Set %d is on %d bucket\n", i, hashes[i][j]);
      // }
    }
  }

  printf("\n\n=== Buckets on each stage ===\n");
  printMatrix(stages, buckets, counts);

  deallocateMatrix(stages, counts);
}

typedef struct initialDataType {
  int nSets;
  int setSize;
  int stages;
  int buckets;
} initialDataType;

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
    initialData->stages = 2;
    initialData->buckets = 4;

    // Dataset params ///
    initialData->nSets = 10;
    initialData->setSize = 4;
    
    for(i = 0; i < nProcess; i++) {
      MPI_Send(initialData, 1, mpiInitialData, i, 0, MPI_COMM_WORLD);
    }
  } else {
    MPI_Recv(initialData, 1, mpiInitialData, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // printf(
    //   "Proccess[%d] receiving stages %d buckets %d nSets %d setSize %d\n", 
    //   myRank,
    //   initialData->stages,
    //   initialData->buckets,
    //   initialData->nSets,
    //   initialData->setSize
    // );
  }
}

//TODO: FORGOT THAT IT NEEDS TO RUN IN DIFFERENT MACHINES
int main () {
  // Parallel variables //
  int nProcess, myRank, i, stages, buckets;
  double start, end;

  // Initialize the MPI enviroment //
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcess);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  initialDataType initialData;
  MPI_Datatype mpiInitialData = createMpiInitialStruct();
  setInitialDataOnProccess(&initialData, &start, mpiInitialData, myRank, nProcess);

  // Generating dataset //
  
  // Size of each partition //
  int partitionSize = nProcess > 0 ? ceil(initialData.nSets/nProcess) : initialData.nSets;

  int *sets = allocateVector(partitionSize * initialData.setSize);
  if(myRank == 0) {
    int *allSets = allocateVector(initialData.nSets * initialData.setSize);
    generateRandomSets(initialData.nSets, initialData.setSize, allSets);

    memcpy(sets, allSets, partitionSize * initialData.setSize * sizeof(int));

    // printf("ALL SETS\n");
    // printVectorAsMatrix(initialData.nSets, initialData.setSize, allSets);

    for(i = 1; i < nProcess; i++) { 
      MPI_Send(&allSets[i * partitionSize * initialData.setSize], partitionSize * initialData.setSize, MPI_INT, i, 0, MPI_COMM_WORLD);
    }

    free(allSets);
  } else {
    MPI_Recv(sets, partitionSize * initialData.setSize, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // printf("SET PROCCESS\n");
    // printVectorAsMatrix(partitionSize, initialData.setSize, sets);
    // printf("\n");
  }

  // Compute hash function coefficients //
  int signatureSize = getSignatureSize(initialData.stages);
  
  // Generating hash coefs //
  int *coefs = allocateVector(signatureSize*2);
  if(myRank == 0) {
    for(i = 0; i < signatureSize * 2; i++) {
      coefs[i] = (rand() % LARGE_PRIME) + 1;
    }

    // printf("==== Coefs ====\n");
    // printVectorAsMatrix(signatureSize, 2, coefs);
    for(i = 0; i < nProcess; i++) {
      MPI_Send(coefs, signatureSize*2, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
  } else {
    MPI_Recv(coefs, signatureSize*2, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  
  // printf("Rank[%d]: Starts on %d - Ends on %d\n", myRank, partitionStart, partitionEnd);

  // Generating hashes //
  int hashesSize = partitionSize * initialData.stages;
  if(myRank == 0) {
    // Spliting sets between proccess //
    int partitionStart = myRank * partitionSize;
    int partitionEnd = partitionStart + partitionSize;
    int *hashes = allocateVector(initialData.nSets * initialData.stages);
    
    hashDataset(
      hashes,
      partitionSize,
      initialData.setSize,
      sets,
      initialData.stages,
      initialData.buckets,
      coefs,
      signatureSize,
      partitionStart,
      partitionEnd
    );

    for(i = 1; i < nProcess; i++) {
      printf("Putting on position %d\n", i * partitionSize);
      int *partitionHashes = allocateVector(hashesSize);
      MPI_Recv(partitionHashes, hashesSize, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printVectorAsMatrix(partitionSize, initialData.stages, partitionHashes);
      memcpy(&hashes[hashesSize * i], partitionHashes, hashesSize * sizeof(int));
      // printVectorAsMatrix(initialData.nSets, initialData.stages, hashes);
    }

    printf("\n=====HASHES=====\n");
    printVectorAsMatrix(initialData.nSets, initialData.stages, hashes);

    printElementsPerBucket(hashes, initialData.nSets, initialData.stages, initialData.buckets);
  } else {
    int *hashes = allocateVector(hashesSize);
    // printVectorAsMatrix(partitionSize, initialData.setSize, sets);
    hashDataset(
      hashes,
      partitionSize,
      initialData.setSize,
      sets,
      initialData.stages,
      initialData.buckets,
      coefs,
      signatureSize,
      0,
      partitionSize
    );
    
    // printVectorAsMatrix(partitionSize, initialData.stages, hashes);
    MPI_Send(hashes, hashesSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
  }

  // MPI_Barrier(MPI_COMM_WORLD);

  // if(myRank == 0) {
  //   end = MPI_Wtime();
  //   printf("%.6lf\n", (end - start)*1000.0);
  // }

  MPI_Finalize();

  // deallocateMatrix(nSets, hashes);
  // deallocateMatrix(signatureSize, coefs);
  // deallocateMatrix(nSets, sets);
	return 0;
}