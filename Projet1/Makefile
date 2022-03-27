#MACROS
CC = gcc
CFLAGS = -Wall -g -lm -fopenmp
TARGET = radixsort

all: ${TARGET}

radixsort.o: radixsort.c
	${CC} -c -o radixsort.o radixsort.c ${CFLAGS}

${TARGET} : radixsort.o
	${CC} -o ${TARGET} radixsort.o ${CFLAGS}

clean:
	-rm -f ${TARGET} *~ *.o *.bak core \#*