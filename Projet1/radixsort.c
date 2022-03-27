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

/* This method allocate a struct array which is a struct containing a tab of int and it's size
*/
struct array * allocateArray(int size) {
    struct array * tmp = malloc(sizeof(struct array));
    tmp->size = size;
    tmp->tab = malloc(size*sizeof(int));
    return tmp;
}

/*This method copies the elements of an array source inside the array dest. 
* Important : source and dest must have the same size !
*/
void copy_array(struct array* source, struct array* dest) {
    #pragma omp parallel for
	for(int i=0; i<source->size; i++) {
		dest->tab[i] = source->tab[i];
	}	
}

/*This method copies the elements of an array source at the end of the array dest. 
* Important : dest must have a size which is twice the size of source !
*/
void copy_array_double(struct array* source, struct array* dest) {
    #pragma omp parallel for
	for(int i=0; i<source->size; i++) {
		dest->tab[i+source->size] = source->tab[i];
	}	
}

/* This method copies the elements of an array source inside of the array dest 
*  Important : start and end must be reachable indices inside source and the substrat of end and start must be the size of dest !
*/
void copy_array_part(struct array* source, struct array* dest, int start, int end) {
    #pragma omp parallel for
    for(int i=start; i<end; i++) {
		dest->tab[i-start] = source->tab[i];
	}	
}

/* This method prints the content of an array and its size on the standard output. Useful for debug !
*/
void printArray(struct array * tmp) {
	printf("---- Array of size %i ---- \n", tmp->size);
	int size = tmp->size;
	int i;
	for (i = 0; i < size; ++i) {
		printf("%i ", tmp->tab[i]);
	}
	printf("\n");
}

/* This method prints the content of an array on the standard output. 
*  Used to printed the content of the initial array and the sorted array inside the output file.
*/
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

/* This method fills an array array of size arraySize with random integers between O and N. 
*  Used when there is no input file specified on the command line.
*/
void fillArray(int N, int arraySize, struct array *array) {
    for (int i=0; i<array->size; i++) {
        array->tab[i] = rand()%(N+1);
    }
}

/* This method fills an array array of size arraySize with the integers inside the input file.
*  This method can produce errors is the content of the file is not suitable with the arguments of the command line
*  (bad input file, to much or to few integers, integers larger than N, ...)
*  Used when there is an input file specified on the command line.
*/
void fillArrayFromFile(int N, int arraySize, char *inputFile, struct array *array) {
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

/* This method formats the content that will be print inside the output file. I redirect the standard output to the file print inside
* it and redirect the standart output to the terminal. 
*/
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

/* This method translates an integer in base 10 in binary and store the bits in a reverse ordre inside the array binary
*  Example : decToBinary(10,binary[]) -> binary = {0,1,0,1}
*/
void decToBinary(int n, struct array *binary) {
    int i=0;
    while(n>0) {
        binary->tab[i] = n%2;
        n = n/2;
        i++;
    }
}

/* This method compute the ascent method for both prefix and suffix on an array which is already size*2 of original array and have
*  integers at the end.
*  This method adds sub elements and store the result inside the node
*  No problem with parallelisation, the values of the nodes didn't change
*  Used in prefix() and suffix()
*  Important : original must have a size equals to 2^n
*  Example : original = {0,0,0,0,1,2,3,4} => original = {0,10,3,7,1,2,3,4}
*/
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

/* This method compute the descent method for prefix on an array which is already size*2 of original array and have
*  been ascent.
*  This method adds the parent and left child and store the result inside the right child, it then stores the value
*  of the parent inside the left child
*  Problem with parallelisation, see report
*  Used in prefix()
*  Important : original must have a size equals to 2^n
*  Example : original = {0,10,3,7,1,2,3,4} => original = {0,0,0,3,0,1,3,6}
*/
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

/* This method compute the descent method for suffix on an array which is already size*2 of original array and have
*  been ascent.
*  This method adds the parent and right child and store the result inside the left child, it then stores the value
*  of the parent inside the right child
*  Problem with parallelisation, see report
*  Used in suffix()
*  Important : original must have a size equals to 2^n
*  Example : original = {0,10,3,7,1,2,3,4} => original = {0,0,7,0,9,7,4,0}
*/
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

/* This method applies the prefix method on the source array and store the result inside the array result
*  It first copies the element of source at the end of result using copy_array_double(). It then computes ascent() of the array result 
*  It finally computes the prefix_descent() on it.
*  Used in prefix_format()
*  Important : result must have a size which is twice the size of source and both must be equals to 2^n !
*  Example : source = {1,2,3,4} and result = {0,0,0,0,0,0,0,0} => source = {1,2,3,4} and result = {0,0,0,3,0,1,3,6}
*/ 
void prefix(struct array *source,struct array *result) {
    copy_array_double(source, result);
    ascent(result);
    prefix_descent(result);
} 

/* This method applies the suffix method on the source array and store the result inside the array result
*  It first copies the element of source at the end of result using copy_array_double(). It then computes ascent() of the array result 
*  It finally computes the suffix_descent() on it.
*  Used in suffix_format()
*  Important : result must have a size which is twice the size of source and both must be equals to 2^n !
*  Example : source = {1,2,3,4} and result = {0,0,0,0,0,0,0,0} => source = {1,2,3,4} and result = {0,0,7,0,9,7,4,0}
*/ 
void suffix(struct array *source, struct array *result) {
    copy_array_double(source, result);
    ascent(result);
    suffix_descent(result);
}

/* This method compute the prefix method on source and store the result inside result.
*  First, it creates a temporary array of size twice the source size, then it computes the prefix operation on it, it copies
*  the end of the temporary array inside the result array and finally releases the memory allocated for this operation.
*  Important : source and result must have the same size and both must be equals to 2^n !
*  Example : source = {1,2,3,4} and result = {0,0,0,0} => source = {1,2,3,4} and result = {0,1,3,6}
*/
void prefix_format(struct array* source, struct array* result) {
    struct array* tmp = allocateArray(source->size*2);
    prefix(source,tmp);
    #pragma omp parallel for
    for (int i=0; i<source->size; i++) {
        result->tab[i] = tmp->tab[i+source->size];
    }
    free(tmp->tab);
    free(tmp);
}

/* This method compute the suffix method on source and store the result inside result.
*  First, it creates a temporary array of size twice the source size, then it computes the suffix operation on it, it copies
*  the end of the temporary array inside the result array and finally releases the memory allocated for this operation.
*  Important : source and result must have the same size and both must be equals to 2^n !
*  Example : source = {1,2,3,4} and result = {0,0,0,0} => source = {1,2,3,4} and result = {9,7,4,0}
*/
void suffix_format(struct array* source, struct array* result) {
    struct array* tmp = allocateArray(source->size*2);
    suffix(source,tmp);
    #pragma omp parallel for
    for (int i=0; i<source->size; i++) {
        result->tab[i] = tmp->tab[i+source->size];
    }
    free(tmp->tab);
    free(tmp);
}

/* This method compute the prefix method on source and store the result inside result but with a size which is not 2^n.
*  The idea is to split the source array into sub array of size 2^n starting from the larger size possible, applies the
*  prefix operation on it and copies the result inside the result array.
*  First, it creates a temporary array of size the number of bits needed to write the size of source in binary, then it 
*  converts the source size into bits and stores the result inside bits. Then, it scans the bits array, if the bit is 1, 
*  it creates the source and result subarrays, copies the corresponding source part inside source subarray using copy_array_part,
*  and applies prefix_format() on it. 
*  If finally copies the elements of the sub array result at the corresponding index inside the result array adding max and 
*  prevsource. Max is the latest element of the previous result subarray and prevsource is the latest element of the source subarray.
*  It sets max and prevsource, releases the memory allocated and computes another loop. At the very end, it releases the memory 
*  allocated for the bits.
*  Used in splitGen()
*  Important : source and result must have the same size !
*  Example : source = {1,2,3,4,5,6} and result = {0,0,0,0,0,0} => source = {1,2,3,4,5,6} and result = {0,1,3,6,10,15}
*/
void prefixGen(struct array *source, struct array *result) {
    int size = source->size; 
    int bitSize = (int) (log2(size)) +1 ;
    struct array* bits = allocateArray(bitSize);
    decToBinary(size,bits);
    int max = 0;
    int start = 0;
    int end;
    int prevsource = 0;
    for (int i=bits->size-1; i>=0; i--) {
        if (bits->tab[i]) {
            end = start + (1<<i);
            struct array *sourcetmp = allocateArray(1<<i);
            struct array *resulttmp = allocateArray(1<<i);
            copy_array_part(source,sourcetmp,start,end);
            prefix_format(sourcetmp,resulttmp);
            for(int j=0; j<resulttmp->size; j++) {
                result->tab[j+start] = resulttmp->tab[j] + max + prevsource;
            }
            max = result->tab[end-1];
            prevsource = source->tab[end-1];
            free(sourcetmp->tab);
            free(resulttmp->tab);
            free(sourcetmp);
            free(resulttmp);
        start = end;
        }
    }
    free(bits->tab);
    free(bits);
}

/* This method compute the suffix method on source and store the result inside result but with a size which is not 2^n.
*  The idea is to split the source array into sub array of size 2^n starting from the larger size possible, applies the
*  suffix operation on it and copies the result inside the result array.
*  First, it creates a temporary array of size the number of bits needed to write the size of source in binary, then it 
*  converts the source size into bits and stores the result inside bits. Then, it scans the bits array, if the bit is 1, 
*  it creates the source and result subarrays, copies the corresponding source part inside source subarray using copy_array_part,
*  and applies suffix_format() on it. 
*  If finally copies the elements of the sub array result at the corresponding index inside the result array adding max and 
*  prevsource. Max is the first element of the previous result subarray and prevsource is the first element of the source subarray.
*  It sets max and prevsource, releases the memory allocated and computes another loop. At the very end, it releases the memory 
*  allocated for the bits.
*  For the suffixGen, it is important to add the source element to the result element in order to delete the gap for the down array.
*  Used in splitGen()
*  Important : source and result must have the same size !
*  Example : source = {1,2,3,4,5,6} and result = {0,0,0,0,0,0} => source = {1,2,3,4,5,6} and result = {21,20,18,15,11,6}
*/
void suffixGen(struct array *source, struct array *result) {
    int size = source->size; 
    int bitSize = (int) (log2(size)) +1 ;
    struct array* bits = allocateArray(bitSize);
    decToBinary(size,bits);
    int max = 0;
    int start;
    int end = source->size;
    int prevsource = 0;
    for (int i=bits->size-1; i>=0; i--) {
        if (bits->tab[i]) {
            start = end - (1<<i);
            struct array *sourcetmp = allocateArray(1<<i);
            struct array *resulttmp = allocateArray(1<<i);
            copy_array_part(source,sourcetmp,start,end);
            suffix_format(sourcetmp,resulttmp);
            for(int j=0; j<resulttmp->size; j++) {
                result->tab[j+start] = resulttmp->tab[j] + max + prevsource;
            }
            max = result->tab[start];
            prevsource = sourcetmp->tab[0];
            free(sourcetmp->tab);
            free(resulttmp->tab);
            free(sourcetmp);
            free(resulttmp);
        end = start;
        }
    }
    free(bits->tab);
    free(bits);
    #pragma omp parallel for
    for (int i=0; i<source->size; i++) {
        result->tab[i] += source->tab[i];
    }
} 

/* This method creates the up array used in splitGen(), it just subtract the size of the array from each element of the source
*  No problem with parallelisation, each operation is unique.
*/
void createUp(struct array* source) {
    #pragma omp parallel for 
    for(int i=0; i<source->size; i++) {
        source->tab[i] = source->size - source->tab[i];
    }
}

/* This method reverse the bits inside the flags array.
*  Used in splitGen()
*/
void notFlags(struct array* flags) {
    #pragma omp parallel for
    for (int i=0; i<flags->size; i++) {
        flags->tab[i] = !flags->tab[i];
    }
}

/* This method stores the elements of the source array inside result array following the indices inside index array.
*  Important : source, index and result arrays must have the same size and each element of index must be unique !
*/
void permute(struct array *source, struct array *index, struct array *result) {
    #pragma omp parallel for
    for (int i=0; i<source->size; i++) {
        int indice = index->tab[i];
        result->tab[indice] = source->tab[i]; 
    }
}

/* This method applies the split methods on source array and stores the result inside result array following the flags.
*  It allocates memory for the index, up and down array, computes the operation suffix, and prefix on flags and not flags.
*  It then deduces the values of the index elements following the bit inside flags, realizes the permute operation and
*  finally release the memory allocated.
*  Important : source, flags and result arrays must have the same size !
*/
void splitGen(struct array *source, struct array *flags, struct array *result) {
    struct array *index = allocateArray(source->size);
    struct array *up = allocateArray(source->size);
    struct array *down = allocateArray(source->size);     
    
    suffixGen(flags,up);
    createUp(up);
    notFlags(flags);
    prefixGen(flags,down);
    notFlags(flags);

    #pragma omp parallel for 
    for (int i=0; i<source->size; i++) {
        if (flags->tab[i]) {
            index->tab[i] = up->tab[i];
        }
        else {
            index->tab[i] = down->tab[i];
        }
    }
    permute(source,index,result);

    free(index->tab);
    free(up->tab);
    free(down->tab);
    free(index);
    free(up);
    free(down);
}

/* This method retrieve the indice's bit of each elements of the source array and stores the result inside flags
*  Important : source and flags must have the same size and indice must be a value between 0 and source size-1.
*/
void generateFlags(struct array *source, struct array *flags, int indice) {
    #pragma omp parallel for 
    for(int i=0; i<source->size; i++) {
        flags->tab[i] = ((source->tab[i] & (1<<indice)) != 0);
    }
}

/* This method is the one that perform the radix sort.
*  It allocates memory for the flags and tmp arrays. It then loops the number of times it is necessary to write N in bits. 
*  Tmp is an array used to stores the values of the result of splitGen(). 
*  It finally release the memory allocated for flags and tmp array.
*  Important : source and sorted must have the same size !
*/
void radixsort(struct array *source, int N, struct array *sorted) {
    struct array *flags = allocateArray(source->size);
    struct array *tmp = allocateArray(source->size);
    copy_array(source,sorted);
    int i=0;
    while(N>0) {
        generateFlags(sorted,flags,i);
        splitGen(sorted,flags,tmp);
        copy_array(tmp,sorted);
        i++;
        N = N/2;      
    }
    free(flags->tab);
    free(tmp->tab); 
    free(flags);
    free(tmp);
}


/* Main method containing all the variables useful, allocates the memory for the initial and sorted array, check the command line, ...
*  More details in report
*/
int main(int argc, char *argv[]) {
    int N = 16;
    int arraySize = 16;
    char *inputFile = NULL;
    char *outputFile = NULL;
    int nbThreads = 1;
    struct timeval start, end;
    
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
    
    if (N<0 || arraySize<0 || nbThreads<0) {
        printf("Error, please enter positive integers.\nExample : radixsort 1021 12 myfile myoutput 2\n");
        exit(1);
    }

    printf("Launch radix sort with : \nN=%d Size=%d Threads=%d Input=%s Output=%s\n", N, arraySize,nbThreads,inputFile,outputFile);

    struct array *arrayToSort = allocateArray(arraySize);
    struct array *arraySorted = allocateArray(arraySize);

    if (inputFile==NULL) {
        fillArray(N, arraySize, arrayToSort);
    }
    else {
        fillArrayFromFile(N, arraySize, inputFile,arrayToSort);
    }
    
    #ifdef _OPENMP
        omp_set_num_threads(nbThreads);
        omp_set_nested(1); // autorise l'imbrication
        omp_set_max_active_levels(3); // d'un niveau max de 3 par ex
        //printf("Nb th in main:%d - (arbitrarily fixed) maximum nested levels: %d\n", omp_get_num_threads( ),omp_get_max_active_levels());
    #endif

    gettimeofday(&start, NULL);
    radixsort(arrayToSort, N, arraySorted);
	gettimeofday(&end, NULL);
    long executionTime = (end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec);
    printOutput(N,arraySize,nbThreads,executionTime,outputFile, arrayToSort, arraySorted);
    printf("Execution in %ld microseconds\n",executionTime);
    free(arrayToSort->tab);
    free(arraySorted->tab); 
    free(arrayToSort);
    free(arraySorted); 
    printf("Done\n");
}