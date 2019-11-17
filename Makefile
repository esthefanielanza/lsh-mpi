CC = mpicc
# CFLAGS = -g -O0 -Wall -Wextra -Werror -std=c99 -fopenmp
CCFLAGS = -g
LIBS = -lm

all: sequential parallel clean

sequential: sequential.o
	$(CC) $(CFLAGS) lsh_sequential.o -o lsh_sequential $(LIBS)

sequential.o: lsh_sequential.c
	$(CC) $(CFLAGS) -c lsh_sequential.c

parallel: parallel.o
	$(CC) $(CFLAGS) lsh_parallel.o -o lsh_parallel $(LIBS)

parallel.o: lsh_parallel.c
	$(CC) $(CFLAGS) -c lsh_parallel.c

clean:
	rm *.o
