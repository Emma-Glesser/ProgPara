#MACROS
CC = mpicc
CFLAGS = -Wall -g -lm 
TARGET = main_mpi
SOURCE = main

all: ${TARGET}

${SOURCE}.o: ${SOURCE}.c
	${CC} -c -o ${SOURCE}.o ${SOURCE}.c ${CFLAGS}

${TARGET} : ${SOURCE}.o
	${CC} -o ${TARGET} ${SOURCE}.o ${CFLAGS}

clean:
	-rm -f ${TARGET} *~ *.o *.bak core \#*