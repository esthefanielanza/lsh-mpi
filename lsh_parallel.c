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
}

void hashSignature(int *hash, int *signature, int stages, int signatureSize, int buckets) {
  int rows = signatureSize / stages;

  for(int i = 0; i < signatureSize; i++) {
    int stage = min(i / rows, stages - 1);
    hash[stage] = (int)((hash[stage] + (long) signature[i] * LARGE_PRIME) % buckets);
  }
}

void hashDataset(
  int **hashes,
  int nSets,
  int setSize,
  int **sets,
  int stages,
  int buckets,
  int **coefs,
  int signatureSize,
  int partitionStart,
  int partitionEnd
) {
  int i;
  int *set = allocateVector(setSize);
  int *signature = allocateVector(signatureSize);

  // printf("\n\n=== Computing coefficients ===\n");
  // printMatrix(signatureSize, 2, coefs);

  // printf("\n\n=== Calculating Hash ===\n");

  for(i = partitionStart; i < partitionEnd; i++) {
    convertToSet(set, sets[i], setSize);
    // printVector(setSize, set);
    calculateSignature(signature, setSize, signatureSize, set, coefs);
    hashSignature(hashes[i], signature, stages, signatureSize, buckets);
    // printf("Hash[%d]:", i);
    // printf(" : ");
    // printVector(stages, hashes[i]);
  }

  free(set);
  free(signature);
}

void printElementsPerBucket(int **hashes, int nSets, int stages, int buckets) {
  int i, j;
  int **counts = allocateMatrix(stages, buckets);

  // printf("\n\n=== Last stage position ===\n");
  for(i = 0; i < nSets; i++) {
    for(j = 0; j < stages; j++) {
      counts[j][hashes[i][j]]++;
      // if(j == stages - 1) {
      //   printf("Set %d is on %d bucket\n", i, hashes[i][j]);
      // }
    }
  }

  // printf("\n\n=== Buckets on each stage ===\n");
  // printMatrix(stages, buckets, counts);

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
    initialData->setSize = 100;
    
    for(i = 0; i < nProcess; i++) {
      MPI_Send(initialData, 1, mpiInitialData, i, 0, MPI_COMM_WORLD);
    }
  } else {
    MPI_Recv(initialData, 1, mpiInitialData, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf(
      "Proccess[%d] receiving stages %d buckets %d nSets %d setSize %d\n", 
      myRank,
      initialData->stages,
      initialData->buckets,
      initialData->nSets,
      initialData->setSize
    );
  }
}

//TODO: FORGOT THAT IT NEEDS TO RUN IN DIFFERENT MACHINES
int main () {
  // Parallel variables //
  int nProcess, myRank, i, nSets, setSize, stages, buckets;
  double start, end;

  // // Compute hash function coefficients //
  // int signatureSize = getSignatureSize(stages);

  // // Generating dataset //
  // int **sets = allocateMatrix(nSets, setSize);

  // generateRandomSets(nSets, setSize, sets);

  // // printf("=== Sets ===\n");
  // // printMatrix(nSets, setSize, sets);

  // // printf("Signature size %d\n", signatureSize);

  // int **coefs = allocateMatrix(signatureSize, 2);
  // for(i = 0; i < signatureSize; i++) {
  //   coefs[i][0] = (rand() % LARGE_PRIME) + 1;
  //   coefs[i][1] = (rand() % LARGE_PRIME) + 1;
  // }

  // int **hashes = allocateMatrix(nSets, stages);

  // Initialize the MPI enviroment //
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcess);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  initialDataType initialData;
  MPI_Datatype mpiInitialData = createMpiInitialStruct();
  setInitialDataOnProccess(&initialData, &start, mpiInitialData, myRank, nProcess);

  // // Size of each partition //
  // int partitionSize = nProcess > 0 ? ceil(nSets/nProcess) : nSets;

  // // Spliting sets between proccess //
  // int partitionStart = myRank * partitionSize;
  // int partitionEnd = partitionStart + partitionSize;

  // printf("Rank[%d]: Starts on %d - Ends on %d\n", myRank, partitionStart, partitionEnd);

  // // Generating hashes //
  // hashDataset(hashes, nSets, setSize, sets, stages, buckets, coefs, signatureSize, partitionStart, partitionEnd);

  // printElementsPerBucket(hashes, nSets, stages, buckets);

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