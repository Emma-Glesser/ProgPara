#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include <math.h>
#include <ctype.h>
#include <fcntl.h>
#include <unistd.h>


struct array {
	int * tab;
	int size;
};

struct array * allocateArray(int size) {
    struct array * tmp = malloc(sizeof(struct array));
    tmp->size = size;
    tmp->tab = malloc(size*sizeof(int));
    return tmp;
}

void copy_array(struct array* source, struct array* dest) {
	for(int i=0; i<source->size; i++) {
		dest->tab[i] = source->tab[i];
	}	
}

void copy_array_double(struct array* source, struct array* dest) {
	for(int i=0; i<source->size; i++) {
        dest->tab[i] = 0;
		dest->tab[i+source->size] = source->tab[i];
	}	
}

void copy_array_part(struct array* source, struct array* dest, int start, int end) {
    for(int i=start; i<end; i++) {
		dest->tab[i-start] = source->tab[i];
	}	
}

void printArray(struct array * tmp) {
	printf("---- Array of size %i ---- \n", tmp->size);
	int size = tmp->size;
	int i;
	for (i = 0; i < size; ++i) {
		printf("%i ", tmp->tab[i]);
	}
	printf("\n");
}

void printArrayOut(struct array * tmp) {
	int size = tmp->size;
	int i;
	for (i = 0; i < size; ++i) {
			printf("%i ", tmp->tab[i]);
	}
	printf("\n");
}

struct array* fileExampleArray() {
    struct array *result = allocateArray(10);
    result->tab[0] = 5;
    result->tab[1] = 7;
    result->tab[2] = 3;
    result->tab[3] = 1;
    result->tab[4] = 4;
    result->tab[5] = 2;
    result->tab[6] = 7;
    result->tab[7] = 2;
    result->tab[8] = 6;
    result->tab[9] = 0;
    return result;
}

struct array* fileExampleArray2() {
    struct array *result = allocateArray(6);
    result->tab[0] = 3;
    result->tab[1] = 4;
    result->tab[2] = 9;
    result->tab[3] = 7;
    result->tab[4] = 5;
    result->tab[5] = 2;
    return result;
}

void fileArray(int N, int arraySize, struct array *array) {
    for (int i=0; i<array->size; i++) {
        array->tab[i] = rand()%(N+1);
    }
}

void fileArrayFromFile(int N, int arraySize, char *inputFile, struct array *array) {
    int fd = open(inputFile, O_RDONLY);
    int n;
    if (fd ==-1) {
        printf("Error, please enter an existing input file.\n");
        exit(1);
    }
    char c[2]; 
    char prev[1];
    char entry[20]; 
    int nb_entries = 0;
    int index = 0;
    c[0] = '\0'; c[1] = '\0'; prev[0] = ' '; entry[0] = '\0';
    while((n=read(fd,c,1))!=0) {
        if (c[0]>=48 && c[0]<=57) {
            if (prev[0]<48 || prev[0]>57) {
                nb_entries += 1;
                if (nb_entries>arraySize) {
                    printf("Error, to much entries, please enter a correct arraySize.\n");
                    exit(1);
                }
            }
            strcat(entry,c);
            index++;
        }

        if (c[0]<48 || c[0]>57) {
            entry[index] = '\0';
            int nb = atoi(entry);
            if (nb>N) {
                printf("Error, entry larger than N, please enter a correct N.\n");
                exit(1);
            }
            if (index) {
                array->tab[nb_entries-1] = nb;
            }
            entry[0] = '\0';
            index = 0;
        }
        prev[0] = c[0];
    }
    close(fd);
    if(index) {
        int nb = atoi(entry);
        if (nb>N) {
            printf("Error, entry larger than N, please enter a correct N.\n");
            exit(1);
        }
        else {
            array->tab[nb_entries-1] = nb;
        } 
    }
    if (nb_entries<arraySize) {
        printf("Error, to few entries, please enter a correct arraySize.\n");
        exit(1);
    }
}

void printOutput(int N, int size, int nbThreads, long executionTime, char *outputFile, struct array *original, struct array *sorted) {
    int fd = open(outputFile,O_WRONLY|O_CREAT, S_IRWXU|S_IRWXG|S_IRWXO);
    if (fd==-1) {
        printf("Error when opening the outputfile\n");
        exit(1);
    }
    else {
        int old_stdout = dup(1);
        dup2(fd,1);
        printf("N=%d Size=%d Threads=%d Time_in_microsec=%ld\n",N,size,nbThreads,executionTime);
        printArrayOut(original);
        printArrayOut(sorted);
        dup2(old_stdout,1);
    }
}

struct array *decToBinary(int n) {
    int bSize = (int) (log2(n)) +1 ;
    struct array* binary = allocateArray(bSize);
    int i=0;
    while(n>0) {
        binary->tab[i] = n%2;
        n = n/2;
        i++;
    }
    return binary;
}

void ascent(struct array *original) {
    int k = (int) (log2(original->size/2));
    for (int j=(k-1); j>=0; j--) {
        int start = 1<<j;
        int end = (1<<(j+1)) - 1;
        #pragma omp parallel for
        for (int i=start; i<=end; i++) {
            original->tab[i] = original->tab[2*i] + original->tab[2*i+1];
        }
    }
}

void prefix_descent(struct array *original) {
    original->tab[1]=0;
	int k = (int) (log2(original->size/2));
    for(int j=0; j<k; j++){
        int start = 1<<j;
        int end = 1<<(j+1);
        #pragma omp parallel for
        for(int i=start; i<end; i++){
            original->tab[2*i+1] = original->tab[i] + original->tab[2*i];
            original->tab[2*i] = original->tab[i];
        }
    }
}

void suffix_descent(struct array *original) {
    original->tab[1]=0;
	int k = (int) (log2(original->size/2));
    for(int j=0; j<k; j++){
        int start = 1<<j;
        int end = 1<<(j+1);
        #pragma omp parallel for
        for(int i=start; i<end; i++){
            original->tab[2*i] = original->tab[i] + original->tab[2*i+1];
            original->tab[2*i+1] = original->tab[i];

        }
    }
}

void prefix(struct array *source,struct array *result) {
    copy_array_double(source, result);
    ascent(result);
    prefix_descent(result);
} 

void suffix(struct array *source, struct array *result) {
    copy_array_double(source, result);
    ascent(result);
    suffix_descent(result);
}

struct array* scan(struct array *source) {
    struct array* tmp = allocateArray(source->size*2);
    struct array* result = allocateArray(source->size);
    prefix(source, tmp);
    #pragma omp parallel for
    for (int i=0; i<source->size; i++) {
        result->tab[i] = tmp->tab[i+source->size];
    }
    //free tmp
    return result;
}

void prefix_format(struct array* source, struct array* result) {
    struct array* tmp = allocateArray(source->size*2);
    prefix(source,tmp);
    #pragma omp parallel for
    for (int i=0; i<source->size; i++) {
        result->tab[i] = tmp->tab[i+source->size];
    }
    //free tmp
}

void suffix_format(struct array* source, struct array* result) {
    struct array* tmp = allocateArray(source->size*2);
    suffix(source,tmp);
    #pragma omp parallel for
    for (int i=0; i<source->size; i++) {
        result->tab[i] = tmp->tab[i+source->size];
    }
    //free tmp
}

void prefixGen(struct array *source, struct array *result) {
    int size = source->size; 
    struct array *sizeBinary = decToBinary(size);
    int max = 0;
    int start = 0;
    int end;
    int prevsource = 0;
    for (int i=sizeBinary->size-1; i>=0; i--) {
        if (sizeBinary->tab[i]) {
            end = start + (1<<i);
            struct array *sourcetmp = allocateArray(1<<i);
            struct array *resulttmp = allocateArray(1<<i);
            copy_array_part(source,sourcetmp,start,end);
            prefix_format(sourcetmp,resulttmp);
            for(int j=0; j<resulttmp->size; j++) {
                result->tab[j+start] = resulttmp->tab[j] + max + prevsource;
            }
            max = resulttmp->tab[(1<<i)-1];
            prevsource = source->tab[(1<<i)-1];
        start = end;
        }
    }
    //free tmpsource et tmpresult
}

void suffixGen(struct array *source, struct array *result) {
    int size = source->size; 
    struct array *sizeBinary = decToBinary(size);
    int max = 0;
    int start;
    int end = source->size;
    int prevsource = 0;
    for (int i=sizeBinary->size-1; i>=0; i--) {
        if (sizeBinary->tab[i]) {
            start = end - (1<<i);
            struct array *sourcetmp = allocateArray(1<<i);
            struct array *resulttmp = allocateArray(1<<i);
            copy_array_part(source,sourcetmp,start,end);
            suffix_format(sourcetmp,resulttmp);
            for(int j=0; j<resulttmp->size; j++) {
                result->tab[j+start] = resulttmp->tab[j] + max + prevsource;
            }
            max = resulttmp->tab[0];
            prevsource = sourcetmp->tab[0];
        end = start;
        }
    }
    //free tmpsource et tmpresult
    #pragma omp parallel for
    for (int i=0; i<source->size; i++) {
        result->tab[i] += source->tab[i];
    }
} 

struct array* createUpArray(struct array *source) {
    struct array *tmp = allocateArray(source->size*2);
    struct array *result = allocateArray(source->size);
    suffix(source,tmp);
    #pragma omp parallel for
    for(int i=0; i<source->size; i++) {
        result->tab[i] = source->size - (tmp->tab[i+source->size] + source->tab[i]);
    }
    //free tmp
    return result;
}

void createUp(struct array* source) {
    for(int i=0; i<source->size; i++) {
        source->tab[i] = source->size - source->tab[i];
    }
}

struct array* not(struct array* flags) {
    struct array *result = allocateArray(flags->size);
    #pragma omp parallel for
    for (int i=0; i<flags->size; i++) {
        result->tab[i] = !flags->tab[i];
    }
    return result;
}

struct array* permute(struct array *source, struct array *index) {
    struct array *result = allocateArray(source->size);
    #pragma omp parallel for
    for (int i=0; i<source->size; i++) {
        result->tab[index->tab[i]] = source->tab[i]; 
    }
    return result;
}

struct array* split(struct array *source, struct array *flags) {
    struct array *index = allocateArray(source->size);
    struct array *up = createUpArray(flags);
    struct array *down = scan(not(flags)); 

    #pragma omp parallel for 
    for (int i=0; i<source->size; i++) {
        if (flags->tab[i]) {
            index->tab[i] = up->tab[i];
        }
        else {
            index->tab[i] = down->tab[i];
        }
    }
    struct array* result = permute(source,index);
    return result;
}

struct array* splitGen(struct array *source, struct array *flags) {
    struct array *index = allocateArray(source->size);
    struct array *up = allocateArray(source->size);
    struct array *down = allocateArray(source->size); 
    suffixGen(flags,up);
    createUp(up);
    prefixGen(not(flags),down);

    #pragma omp parallel for 
    for (int i=0; i<source->size; i++) {
        if (flags->tab[i]) {
            index->tab[i] = up->tab[i];
        }
        else {
            index->tab[i] = down->tab[i];
        }
    }
    struct array* result = permute(source,index);
    return result;
}

struct array* generateFlags(struct array *source, int indice) {
    struct array *flags = allocateArray(source->size);
    for(int i=0; i<source->size; i++) {
        flags->tab[i] = ((source->tab[i] & (1<<indice)) != 0);
    }
    return flags;
}

struct array* radixsort(struct array *source, int N) {
    struct array *sorted = allocateArray(source->size);
    copy_array(source,sorted);
    int i=0;
    while(N>0) {
        struct array *flags = generateFlags(sorted,i);
        sorted = splitGen(sorted,flags);
        i++;
        N = N/2;
    }
    return sorted;
}

int main(int argc, char *argv[]) {
    int N = 16;
    int arraySize = 16;
    char *inputFile = NULL;
    char *outputFile = NULL;
    int nbThreads = 1;
    //struct timeval start, end;
    
    if (argc<4 || argc>6) {
        printf("Error, missing or to much arguments please enter a correct command.\nExample : radixsort 1021 12 myfile myoutput 2\n");
        exit(1);
    }
    else if (argc==4) {
        N = atoi(argv[1]);
        arraySize = atoi(argv[2]);
        outputFile = argv[3];
    }
    else if (argc==6) {
        N = atoi(argv[1]);
        arraySize = atoi(argv[2]);
        inputFile = argv[3];
        outputFile = argv[4];
        nbThreads = atoi(argv[5]);
    }
    else if (argc==5 && argv[4][0]>=48 && argv[4][0]<=57) {
        N = atoi(argv[1]);
        arraySize = atoi(argv[2]);
        outputFile = argv[3];
        nbThreads = atoi(argv[4]);

    }
    else {
        N = atoi(argv[1]);
        arraySize = atoi(argv[2]);
        inputFile = argv[3];
        outputFile = argv[4];
    }
    
    printf("N=%d Size=%d Threads=%d Input=%s Output=%s\n", N, arraySize,nbThreads,inputFile,outputFile);

    struct array *arrayToSort = allocateArray(arraySize);

    if (inputFile==NULL) {
        fileArray(N, arraySize, arrayToSort);
    }
    else {
        fileArrayFromFile(N, arraySize, inputFile,arrayToSort);
    }

    struct array *array = fileExampleArray();
    printArray(array);
    struct array *test = radixsort(array,N);
    printArray(test);
    
/*     gettimeofday(&start, NULL);
    struct array *arraySorted = radixsort(arrayToSort, N);
	gettimeofday(&end, NULL);
    long executionTime = (end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec);
    printOutput(N,arraySize,nbThreads,executionTime,outputFile, arrayToSort, arraySorted); */


}