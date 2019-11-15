CC = gcc
CFLAGS = -g -O0 -Wall -Wextra -Werror -std=c99 -fopenmp
LIBS = -lm

all: sequential clean

sequential: sequential.o
	$(CC) $(CFLAGS) lsh_sequential.o -o lsh_sequential $(LIBS)

sequential.o: lsh_sequential.c
	$(CC) $(CFLAGS) -c lsh_sequential.c

# parallel: parallel.o input.o pthreadPool.o
# 	$(CC) $(CFLAGS) parallel.o input.o pthreadPool.o -o parallel $(LIBS)

# parallel.o: parallel.c input.h pthreadPool.h
# 	$(CC) $(CFLAGS) -c parallel.c

# input.o: input.c input.h
# 	$(CC) $(CFLAGS) -c input.c

# pthreadPool.o: pthreadPool.c pthreadPool.h
# 	$(CC) $(CFLAGS) -c pthreadPool.c

clean:
	rm *.o
