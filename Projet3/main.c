#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>

/* LIRE RAPPORT FIN DU CODE LIGNE 489 */

int procID; 
int nbProcess;

/* Structure representant une matrice linéariser par lignes
*/
struct matrix {
    int *array;
    int rows;
    int columns;
};

/* Méthode permettant de récupérer l'élément de la ligne i et de la colonne j de la matrice m 
*/
int get(struct matrix *m, int i, int j) {
    return m->array[j + i*m->columns];
}

/* Méthode permettant de récupérer l'élément i dans le tableau de la matrice m
*/
int getIndex(struct matrix *m, int i) {
    return m->array[i];
}

/* Méthode permettant de changer la valeur de l'élément de la ligne i et de la colonne j de la matrice m
*/
void set(struct matrix *m, int i, int j, int value) {
    m->array[j + i*m->columns]=value;
}

/* Méthode permettant de changer la valeur de l'élément i dans le tableau de la matrice m
*/
void setIndex(struct matrix *m, int i, int value) {
    m->array[i]=value;
}

/* Méthode permettant de récupérer les lignes de la matrice m pour les stocker dans le tableau dest.
 * start et end représentent les bornes des éléments de la matrice m qui doivent être copiés, end n'est pas inclus
 * Important : start doit être supérieur ou égal à 0 et end ne doit pas dépasser le nombre d'éléments de la matrice, le nombre d'élément de 
 *             dest doit être égal à (end-start). 
 * Exemple : matrice = [[1,2,3],[4,5,6],[7,8,9]] donc m = [1,2,3,4,5,6,7,8,9] 
 *           on veut récupérer la ligne 1 et 2 de m donc start = 0, end = 6 et dest = int[6]
 *           getRows(m,dest,0,6) => dest = [1,2,3,4,5,6]
*/
void getRows(struct matrix *m, int *dest, int start, int end) {
    for(int i=0; i<end-start; i++) {
        dest[i] = getIndex(m,start+i);
    }
}

/* Méthode permettant de stocker les éléments du tableau src dans les lignes de la matrice m.
 * start et end représentent les bornes des éléments de la matrice qui vont être changer, end n'est pas inclus.
 * Important : start doit être supérieur ou égal à 0 et end ne doit pas dépasser le nombre d'éléments de la matrice, le nombre d'élément de 
 *             src doit être égal à (end-start). 
 * Exemple : matrice = [[0,0,0],[0,0,0],[0,0,0]] donc m = [0,0,0,0,0,0,0,0,0]
 *           on veut stocker les éléments de src dans les lignes 1 et 2 de m donc start = 0, end = 6 et src = [1,2,3,4,5,6]
 *           setRows(m,src,0,6) => m = [1,2,3,4,5,6,0,0,0]
*/
void setRows(struct matrix *m, int *src, int start, int end) {
    for(int i=0; i<end-start; i++) {
        setIndex(m,start+i,src[i]);
    }
}

/* Méthode qui permet d'allouer l'espace mémoire nécessaire pour une matrice de taille rows * colums. Renvoie une struct matrix *
*/
struct matrix *allocateMatrix(int rows, int columns) {
    struct matrix *m = (struct matrix *) malloc(sizeof(struct matrix));
    m->rows = rows;
    m->columns = columns;
    m->array = (int*) malloc(sizeof(int*) * rows * columns);
    return m;
}

/* Méthode qui permet d'afficher sur stdout la matrice m et ses dimensions
*/
void printMatrix(struct matrix *m) {
    printf("Matrix %d x %d\n",m->rows, m->columns);
    int v;
    for(int i=0; i<m->rows; i++) {
        for (int j=0; j<m->columns; j++) {
            v=get(m,i,j);
            printf("%d ",v);
        }
        printf("\n");
    }
    printf("\n");
}

/* Méthode qui permet d'afficher sur stdout la matrice m mais sans ses dimensions. Utilisée pour écrire le résultat dans le fichier de sortie
*/
void printMatrixOut(struct matrix *m) {
    int v;
    int lenght = m->rows;
    for(int i=0; i<lenght; i++) {
        for (int j=0; j<m->columns; j++) {
            v=get(m,i,j);
            printf("%d ",v);
        }
        if (i!=lenght-1) {
            printf("\n");
        }
    }
}

/* Méthode qui permet de récupérer les dimensions d'une matrice contenue dans le fichier filename et les stocker dans row et column.
 * Utilisée dans main pour récupérer les dimensions de la matrice et du vecteur afin de les allouer pour tous les processus avant de les
 * remplir uniquement dans le processus 0.
*/
void getDimensions(int *row, int *column, char *filename) {
    int fd = open(filename,O_RDONLY);
    if (fd ==-1) {
        printf("Error, please enter an existing input file.\n");
        exit(1);
    }
    char c[1]; 
    int n;
    int rows = 0;
    int columns = 0;
    c[0]='\0';
    while((n=read(fd,c,1))!=0) {
        if (c[0]==' ') {
            columns++;
        }
        if (c[0]=='\n') {
            columns = 0;
            rows++;
        }
    }
    close(fd);
    if (columns!=0) {
        rows++;
    }
    *row = rows;
    *column = columns;
}

/* Méthode qui permet de remplir les éléments de la matrice m avec comme valeur leur indice + 1
*/
void fillMatrix(struct matrix *m) {
    for(int i=0; i<m->rows;i++) {
        for(int j=0; j<m->columns; j++) {
            int value = (j+1) + i*m->columns;
            set(m,i,j,value);
        }
    }
}

/* Méthode qui permet de remplir les éléments de la matrice m avec 0
*/
void fillNullMatrix(struct matrix *m) {
    for(int i=0; i<m->rows;i++) {
        for(int j=0; j<m->columns; j++) {
            set(m,i,j,0);
        }
    }
}

/* Méthode qui permet de remplir les éléments de la matrice m pour construire la matrice identité (que des 0 sauf sur la diagonale)
*/
void fillIDMatrix(struct matrix *m) {
    for(int i=0; i<m->rows;i++) {
        for(int j=0; j<m->columns; j++) {
            if (i==j) {
                set(m,i,j,1);
            }
            else {
                set(m,i,j,0);
            }
        }
    }
}

/* Méthode qui permet de remplir les éléments de la matrice m à partir du fichier filename 
*/
void fillMatrixFromFile(struct matrix *m, char *filename) {
    int fd = open(filename,O_RDONLY);
    if (fd ==-1) {
        printf("Error, please enter an existing input file.\n");
        exit(1);
    }
    char c[2]; 
    char prev[1];
    char entry[20];
    int n;
    int count = 0;
    int index = 0;
    c[0] = '\0'; c[1] = '\0'; prev[0] = ' '; entry[0] = '\0';
    while((n=read(fd,c,1))!=0) {
        if (c[0]>=48 && c[0]<=57) {
            if (prev[0]<48 || prev[0]>57) {
                count++;
            }
            strcat(entry,c);
            index++;
        }
        else {
            entry[index] = '\0'; 
            int nb = atoi(entry);
            if (index) {
                setIndex(m,count-1,nb);
            }
            entry[0] = '\0';
            index = 0;
        }
        prev[0] = c[0];
    }
    close(fd);
}

/* Méthode qui permet d'écrire la matrice result dans le fichier outputfile.
 * Utilisée à la fin pour écrire le résultat du produit matriciel dans le fichier de sortie
*/
void printOutput(char *outputFile, struct matrix *result) {
    int fd = open(outputFile,O_WRONLY|O_CREAT|O_TRUNC, S_IRWXU|S_IRWXG|S_IRWXO);
    if (fd==-1) {
        printf("Error when opening the outputfile\n");
        exit(1);
    }
    else {
        int old_stdout = dup(1);
        dup2(fd,1);
        printMatrixOut(result);
        fflush(stdout);
        dup2(old_stdout,1);
    }
}

/* Méthode qui diffuse le contenu de la matrice m aux autres processus en suivant une boucle. Les processus ne 
 * conservent que la partie qui leur correspond et la stockent dans la matrice dest. 
 * Exemple : Il y a 3 processus et le processus 0 connait la matrice m = [[1,2,3],[4,5,6],[7,8,9]]. Il envoie d'abord
 *           la ligne [7,8,9] au processus 1 qui la reçoit et la transfère au processus 2. Le processus 2 reçoit [7,8,9]
 *           et les stocke dans dest. Le processus 0 envoie ensuite [4,5,6] au 1 qui les stocke dans dest sans renvoyer.
 *           Finalement le processus 0 stocke [1,2,3] dans dest.
 * Important : le nombre d'éléments de dest doit être égale au nombre d'éléments de m divisé par le nombre de processus.
 *             C'est à dire si m est de dimension n*p alors dest doit être de dimensions (n/nbProcess)*p
 * Calcul des coûts de performance voir section scatter dans rapport à la fin ligne 464
*/
void scatter(struct matrix *m, struct matrix *dest) {
    int tag=1;
	MPI_Status status;

    /* Le processus 0 envoie les différents tableaux en commençant par le bas de la matrice puis finit par stocker le haut
     * de la matrice dans dest */
    if(procID==0) {
        int blocksize = (m->rows * m->columns)/nbProcess;
        for(int i=nbProcess-1; i>0; i--) {
            int start = i*blocksize;
            int end = start + blocksize;
            int *block = malloc(sizeof(int)*blocksize);
            getRows(m,block,start,end);
            MPI_Send(block,blocksize,MPI_INT,1%nbProcess,tag,MPI_COMM_WORLD); 
        }
        getRows(m,dest->array,0,blocksize);
    } 

    /* Les autres processus renvoient les messages reçus tant que le nombre de messages reçus est inférieur ou égal à leur rang
     * Ils finissent par stocker leur message destiné dans dest */
    else {
        for(int i=nbProcess-1; i>=procID; i--) {
            MPI_Probe((procID-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
            int number;
            MPI_Get_count(&status, MPI_INT, &number);
            MPI_Recv(dest->array,number,MPI_INT,(procID-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
            if (i != procID) {
                MPI_Send(dest->array,number,MPI_INT,(procID+1)%nbProcess,tag,MPI_COMM_WORLD);
            }
        }
    }
}

/* Méthode qui renvoie le contenu de la matrice src au processus 0 en suivant une boucle. Les processus reçoivent les
 * les messages des processus précédent les renvoie en ajoutant leur propre message. Le processus récupère tous les 
 * messages et stocke leur contenu au bon endroit dans la matrice m.
 * Exemple : Il y a 3 processus et la matrice m est de taille 3*1. Le processus 1 envoie son contenu de src = [2] à 2.
 *           Le processus 2 reçoit le message de 1 et le renvoie à 0 et envoie ensuite son contenu de src = [3]. Le 
 *           processus 0 reçoit ensuite le message contenant [2] et le stocke dans la ligne 2 de m, puis il reçoit [3]
 *           et le stocke dans la ligne 3 de m et finalement stocke son propre résultat [1] dans la première ligne.
 * Important : le nombre d'éléments de src doit être égale au nombre d'éléments de m divisé par le nombre de processus.
 *             C'est à dire si m est de dimension n*p alors m doit être de dimensions (n/nbProcess)*p
 * Calcul des coûts de performance voir section gather dans rapport à la fin ligne 499
*/
void gather(struct matrix *m, struct matrix *src) {
    int tag=10;
	MPI_Status status;

    /* Calcul de la taille du block à recevoir et stocker */
    int blocksize = (m->rows * m->columns)/nbProcess;

    /* Le processus 0 récupère tous les messages en provenance du dernier processus. Il stocke leur contenu dans le tableau
     * temporaire block et puis le stocke dans la matrice grâce à setRows().*/
    if (procID==0) {
        for(int i=1; i<nbProcess; i++) {
            int number;
            int *block = malloc(sizeof(int)*blocksize);
            MPI_Probe((nbProcess-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
            MPI_Get_count(&status, MPI_INT, &number);
            MPI_Recv(block,number,MPI_INT,(nbProcess-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
            int start = i*blocksize ;
            int end = start + blocksize;
            setRows(m,block,start,end);
        }
        setRows(m,src->array,0,blocksize);
    }

    /* L'ordre est respecté car chaque processus attend d'avoir renvoyé tous les messages reçus dans le même ordre avant d'envoyer
     * propre message */
    else {
        for(int i=1; i<procID; i++) {
            int number;
            int *block = malloc(sizeof(int)*blocksize);
            MPI_Probe((procID-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
            MPI_Get_count(&status, MPI_INT, &number);
            MPI_Recv(block,number,MPI_INT,(procID-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
            MPI_Send(block,blocksize,MPI_INT,(procID+1)%nbProcess,tag,MPI_COMM_WORLD); 
        }
        MPI_Send(src->array,blocksize,MPI_INT,(procID+1)%nbProcess,tag,MPI_COMM_WORLD); 
    }
}

/* Méthode qui diffuse le contenu de la matrice m à tous les autres processus en suivant une boucle.
 * Exemple : il y a 3 processus. Le 0 envoie au 1 qui reçoit et envoie au 2 qui reçoit mais n'envoie rien.
 * Calcul des coûts de performance voir section token_ring dans rapport à la fin ligne 536
*/
void token_ring(struct matrix *m) {
	int tag=100;
	MPI_Status status;
    
    /* Le processus 0 envoie le contenu de la matrice m au processus 1 ou lui même si il n'y a qu'un processus*/
    if (procID==0) {
        int size = m->rows * m->columns;
    	MPI_Send(m->array,size,MPI_INT,1%nbProcess,tag,MPI_COMM_WORLD);  
    } 

    /* Les autres processus récupèrent le contenu de la matrice grâce aux processus avant eux, le stockent et le renvoie aux processus 
     * après eux. Ils connaissent la taille du tableau grâce au probe et get_count 
     * Le dernier processus ne renvoie rien.
    */
    else {
		MPI_Probe((procID-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
		int number;
		MPI_Get_count(&status, MPI_INT, &number);
		MPI_Recv(m->array,number,MPI_INT,(procID-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
        if (procID!=nbProcess-1) {
            MPI_Send(m->array,number,MPI_INT,(procID+1)%nbProcess,tag,MPI_COMM_WORLD);
        }
	}
}

/* Méthode qui diffuse les dimensions des matrices A et B sous forme de tableau à tous les autres processus en suivant une boucle.
 * Exemple : il y a 3 processus. Le 0 envoie au 1 qui reçoit et envoie au 2 qui reçoit mais n'envoie rien.
 * Calcul des coûts de performance voir section token_ring_dim dans rapport à la fin ligne 536
*/
void token_ring_dim(int *rA, int *cA, int *rB, int *cB) {
	int tag=2;
	MPI_Status status;
    
    /* Le processus 0 envoie le contenu de la matrice m au processus 1 ou lui même si il n'y a qu'un processus*/
    if (procID==0) {
        int *array = malloc(sizeof(int)*4);
        array[0] = *rA; array[1] = *cA; array[2] = *rB; array[3] = *cB; 
    	MPI_Send(array,4,MPI_INT,1%nbProcess,tag,MPI_COMM_WORLD);  
    } 

    /* Les autres processus récupèrent le contenu de la matrice grâce aux processus avant eux et le renvoie aux processus 
     * après eux. Ils connaissent la taille du tableau grâce au probe et get_count
     * Le dernier processus ne renvoie rien */
    else {
        int number;
        int *array = malloc(sizeof(int)*4);
		MPI_Probe((procID-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
		MPI_Get_count(&status, MPI_INT, &number);
		MPI_Recv(array,number,MPI_INT,(procID-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
        if (procID!=nbProcess-1) {
		    MPI_Send(array,number,MPI_INT,(procID+1)%nbProcess,tag,MPI_COMM_WORLD);
        }
        *rA = array[0]; *cA = array[1]; *rB = array[2]; *cB = array[3]; 
	}
}

/* Méthode qui permet de réaliser un produit matriciel entre la matrixA et la matrixB et qui stocke le résultat dans la matrice result
 * Important : le nombre de colonnes de A doit être égal au nombre de lignes de B et le résultat aura les dimensions (nb de lignes de A) * (nb colonnes de B)
*/
void multiplyMatrix(struct matrix *matrixA, struct matrix *matrixB, struct matrix *result) {
    for(int i=0; i<matrixA->rows; i++) {
        for (int k=0; k<matrixB->columns; k++) {
            int element = 0; 
            for (int j=0; j<matrixA->columns; j++) {
                element += get(matrixA,i,j)*get(matrixB,j,k);
            }
            set(result,i,k,element);
        }
    }
}

/* Méthode qui permet l'appel aux différentes méthodes pour produire le résultat du calcul de la matrice m par le vecteur vector et le stocker
 * dans la matrice result. A ce stade, tous les éléments ont leur place allouée pour les différents processus mais seulement le processus 0
 * a pu remplir les éléments dans la matrice m et dans la matrice vector. 
 * token_ring diffuse le vector aux autres processus
 * scatter diffuse les différentes lignes de la matrice m aux processus correspondants
 * multiplyMatrix éxécute et stocke le résultat du produit matriciel 
 * gather envoie ce résultat au processus 0 qui reconstitue la matrice résultat
*/
void compute(struct matrix *m, struct matrix *vector, struct matrix *result) {
    /* Découpage de la matrice m en sous-matrices de taille (taille m)/nb processus. */
    int submatrixSize = m->rows/nbProcess;
    struct matrix *submatrix = allocateMatrix(submatrixSize,m->columns);
    struct matrix *subresult = allocateMatrix(submatrixSize,vector->columns);
    /* Envoi du vecteur, puis des sous-matrices, calcul du produit matriciel et envoie des résultats*/
    token_ring(vector);
    scatter(m,submatrix);
    multiplyMatrix(submatrix,vector,subresult);
    gather(result,subresult);
}

/* Méthode main qui vérifie les arguments sur la ligne de commande afin de connaitre les fichiers d'entrée pour la matrice et le vecteur
 * Elle calcule ensuite les dimensions des matrices à allouer dans TOUS les processus grâce à la méthode getDimensions, alloue l'espace 
 * mémoire pour ces matrices, et remplie le contenu de la matrice et du vecteur seulement dans le processus 0 grâce à fillMatrixFromFile().
 * Elle appelle la méthode compute pour calculer le produit et écrit ce résultat dans le fichier de sortie "output.txt" uniquement dans le
 * processus 0.
*/
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);
    MPI_Comm_size(MPI_COMM_WORLD, &nbProcess);

    /*Vérification des arguments en ligne de commande et fin de programme si non conformité */
    if(argc<3 || argc>3) {
        if(procID==0) {
            printf("Error there are to much or missing arguments input files for matrix and vector.\nPlease enter a command line as follow : matmult <inputfile matrix> <inputfile vector>\n");
        }
        exit(1);
    }

    /* Initialisation des valeurs */
    char *inputfileMatrix = argv[1];
    char *inputfileVector = argv[2];
    char *outputfile = "./output.txt";
    int rowsM = 0;
    int columnsM = 0;
    int rowsV = 0;
    int columnsV = 0;

    /* On récupère les dimensions des deux matrices dans le processus */
    if (procID==0) {
        getDimensions(&rowsM,&columnsM, inputfileMatrix);
        getDimensions(&rowsV,&columnsV, inputfileVector);
    }

    /* On diffuse aux autres processus */
    token_ring_dim(&rowsM,&columnsM,&rowsV,&columnsV);

    /* On alloue l'espace mémoire dans chaque processus */
    struct matrix *A = allocateMatrix(rowsM,columnsM);
    struct matrix *B = allocateMatrix(rowsV, columnsV);
    struct matrix *C = allocateMatrix(rowsM,columnsV);

    /* On remplit le contenu des matrices seulement s'il s'agit du processus 0 */
    if (procID==0) {
        fillMatrixFromFile(B, inputfileVector);
        fillMatrixFromFile(A, inputfileMatrix);
    }

    /* Produit matriciel */
    compute(A,B,C);

    /* Print des différents matrices et écriture dans le fichier de sortie uniquement pour le processus 0 */
    if (procID==0) {
        printMatrix(A);
        printMatrix(B);
        printMatrix(C);
        printOutput(outputfile,C);
    }

    MPI_Finalize();
    return 0;
}


/* CALCUL DU COUT DES FONCTIONS ET DU COUT TOTAL 

bêta = latence
tau = bande passante 
L = longueur du message

-- Scatter --
    if(procID==0) {
        int blocksize = (m->rows * m->columns)/nbProcess;
        for(int i=nbProcess-1; i>0; i--) {
            int start = i*blocksize;
            int end = start + blocksize;
            int *block = malloc(sizeof(int)*blocksize);
            getRows(m,block,start,end);
            MPI_Send(block,blocksize,MPI_INT,1%nbProcess,tag,MPI_COMM_WORLD); 
        }
        getRows(m,dest->array,0,blocksize);
    } 
    else {
        for(int i=nbProcess-1; i>=procID; i--) {
            MPI_Probe((procID-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
            int number;
            MPI_Get_count(&status, MPI_INT, &number);
            MPI_Recv(dest->array,number,MPI_INT,(procID-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
            if (i != procID) {
                MPI_Send(dest->array,number,MPI_INT,(procID+1)%nbProcess,tag,MPI_COMM_WORLD);
            }
        }
    }

    Le processus 0 envoie (nbProcess-1) messages de longueur L. 
    Les autres processus envoient (nbProcess-1-rang)  messages de longueur L.
    Comme les messages sont envoyés en parallèle finalement c'est le nombre de messages envoyés par le processus 0 qui est retenu
    = max((nbProcess-1),(nbProcess-1-rang)) = (nbProcess-1).
    Chaque message envoyé est de longueur L=blocksize.
    Coût de l'envoi d'un message = (bêta + L*tau).
    Dimension matrice m = n*p
    Finalement coût de la fonction = (nbProcess-1) * (bêta + L*tau)
                                   = (nbProcess-1) * (bêta + blocksize*tau)
                                   = (nbProcess-1) * (bêta + (n*p/nbProcess)*tau)

-- Gather -- 
    int blocksize = (m->rows * m->columns)/nbProcess;
    if (procID==0) {
        for(int i=1; i<nbProcess; i++) {
            int number;
            int *block = malloc(sizeof(int)*blocksize);
            MPI_Probe((nbProcess-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
            MPI_Get_count(&status, MPI_INT, &number);
            MPI_Recv(block,number,MPI_INT,(nbProcess-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
            int start = i*blocksize ;
            int end = start + blocksize;
            setRows(m,block,start,end);
        }
        setRows(m,src->array,0,blocksize);
    }
    else {
        for(int i=1; i<procID; i++) {
            int number;
            int *block = malloc(sizeof(int)*blocksize);
            MPI_Probe((procID-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
            MPI_Get_count(&status, MPI_INT, &number);
            MPI_Recv(block,number,MPI_INT,(procID-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
            MPI_Send(block,blocksize,MPI_INT,(procID+1)%nbProcess,tag,MPI_COMM_WORLD); 
        }
        MPI_Send(src->array,blocksize,MPI_INT,(procID+1)%nbProcess,tag,MPI_COMM_WORLD); 
    }

    Le processus 0 n'envoie aucun message.
    Les autres processus envoient chaque un message plus tous ceux qu'ils reçoivent donc chaque processus autre que 0 envoie (rang) messages.
    Comme les messages sont envoyés en parallèle, c'est le nombre de messages envoyés par le dernier processus qui compte = max(rang) = (nbProcess-1).
    Chaque message envoyé est de longueur L=blocksize.
    Coût de l'envoi d'un message = (bêta + L*tau)
    Dimension matrice m = n*p
    Finalement coût de la fonction = (nbProcess-1) * (bêta + L*tau)
                                   = (nbProcess-1) * (bêta + blocksize*tau)
                                   = (nbProcess-1) * (bêta + (n*p/nbProcess)*tau)

-- Token ring -- 
    if (procID==0) {
        int size = m->rows * m->columns;
    	MPI_Send(m->array,size,MPI_INT,1%nbProcess,tag,MPI_COMM_WORLD);  
    } 
    else {
		MPI_Probe((procID-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
		int number;
		MPI_Get_count(&status, MPI_INT, &number);
		MPI_Recv(m->array,number,MPI_INT,(procID-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
        if (procID!=nbProcess-1) {
            MPI_Send(m->array,number,MPI_INT,(procID+1)%nbProcess,tag,MPI_COMM_WORLD);
        }
	}

    Le processus 0 envoie un message de longueur L au processus 1.
    Le dernier processus n'envoie pas de message.
    Les autres processus envoient un message chacun de longueur L au processus suivant.
    Finalement il y aura (nbProcess-1) messages envoyés mais comme ils sont tous envoyés en parallèle, on ne compte que le coût d'un seul message.
    Chaque message envoyé est de longueur L=taille_matrice.
    Coût de l'envoi d'un message = (bêta + L*tau)
    Dimension matrice m = n*p
    Finalement coût de la fonction = 1 * (bêta + L*tau)
                                   = (bêta + blocksize*tau)
                                   = (bêta + n*p*tau)

-- Token Ring Dim --
    if (procID==0) {
        int *array = malloc(sizeof(int)*4);
        array[0] = *rA; array[1] = *cA; array[2] = *rB; array[3] = *cB; 
    	MPI_Send(array,4,MPI_INT,1%nbProcess,tag,MPI_COMM_WORLD);  
    } 
    else {
        int number;
        int *array = malloc(sizeof(int)*4);
		MPI_Probe((procID-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
		MPI_Get_count(&status, MPI_INT, &number);
		MPI_Recv(array,number,MPI_INT,(procID-1)%nbProcess,tag,MPI_COMM_WORLD,&status);
        if (procID!=nbProcess-1) {
		    MPI_Send(array,number,MPI_INT,(procID+1)%nbProcess,tag,MPI_COMM_WORLD);
        }
        *rA = array[0]; *cA = array[1]; *rB = array[2]; *cB = array[3]; 
	}

    Le processus 0 envoie un message de longueur L au processus 1.
    Le dernier processus n'envoie pas de message.
    Les autres processus envoient un message chacun de longueur L au processus suivant.
    Finalement il y aura (nbProcess-1) messages envoyés mais comme ils sont tous envoyés en parallèle, on ne compte que le coût d'un seul message.
    Chaque message envoyé est de longueur L=4 car on envoie que 4 entiers, les dimensions des 2 matrices.
    Coût de l'envoi d'un message = (bêta + L*tau)
    Finalement coût de la fonction = 1 * (bêta + L*tau)
                                   = (bêta + 4*tau)
                                   = (bêta + 4*tau)

-- Coût total du programme --

    Matrice m du programme de dimension n * p
    Vecteur vector du programme de dimension p * q
    Résultat resutl du programme de dimesion n * q

    int submatrixSize = m->rows/nbProcess;
    struct matrix *submatrix = allocateMatrix(submatrixSize,m->columns);
    struct matrix *subresult = allocateMatrix(submatrixSize,vector->columns);
    token_ring(vector);
    scatter(m,submatrix);
    multiplyMatrix(submatrix,vector,subresult);
    gather(result,subresult);

    Coût calculé scatter = (nbProcess-1) * (bêta + (n*p/nbProcess)*tau) 
    On remplace par les dimensions de m => (nbProcess-1) * (bêta + (n*p/nbProcess)*tau) 

    Coût calculé gather = (nbProcess-1) * (bêta + (n*p/nbProcess)*tau)
    On remplace par les dimensions de result => (nbProcess-1) * (bêta + (n*q/nbProcess)*tau)

    Coût calculé token_ring = (bêta + n*p*tau)
    On remplace par les dimensions de vector => (bêta + p*q*tau)

    Coût calculé token_ring_dim = (bêta + 4*tau)

    Coût total = scatter + gather + token_ring + token_ring_dim
               = (nbProcess-1) * (bêta + (n*p/nbProcess)*tau) + (nbProcess-1) * (bêta + (n*q/nbProcess)*tau) + (bêta + p*q*tau) + (bêta + 4*tau)
               = (nbProcess-1) * (bêta + (n*p/nbProcess)*tau + bêta + (n*q/nbProcess)*tau) + (bêta + p*q*tau) + (bêta + 4*tau)
               = (nbProcess-1) * (2*bêta + (n*p/nbProcess)*tau + (n*q/nbProcess)*tau) + (bêta + p*q*tau) + (bêta + 4*tau)
               = (nbProcess-1) * (2*bêta + tau*(n*p/nbProcess + n*q/nbProcess)) + (bêta + p*q*tau) + (bêta + 4*tau)
               = (nbProcess-1) * (2*bêta + tau*n*(p/nbProcess + q/nbProcess)) + (bêta + p*q*tau) + (bêta + 4*tau)
               = (nbProcess-1) * (2*bêta + tau*n*(p+q)/nbProcess)) + (bêta + p*q*tau) + (bêta + 4*tau)
               = (nbProcess-1) * (2*bêta + tau*n*(p+q)/nbProcess)) + (2*bêta + (p*q+4)*tau)

    Dans notre cas q = 1 donc :
    Coût total = (nbProcess-1) * (2*bêta + tau*n*(p+1)/nbProcess)) + (2*bêta + (p+4)*tau)

*/