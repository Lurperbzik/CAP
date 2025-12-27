#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h> 

#define R 0.1
#define ESTACIONARIO 1e-6
#define LOG_FILE "resultados.csv"
#define INF 99999

void* hilo_func(void *);

struct Matriz {
    unsigned int nfil, ncol;
    float **actual;
    float **anterior;
};
struct Matriz matriz = {0};

float dt = INF;
int nhilos;
struct Status {
    double tiempo;
    long iteraciones_temporales;
};
struct Status status = {0};

int main(int argc, char * argv[]) {
    /*
        VARIABLES LOCALES
    */
    int i;
    long id;
    struct timeval start, end;
    FILE *f;
    

    /*
        LECTURA DE ARGUMENTOS
    */
    if (argc < 4) {
        printf("Error de uso: %s [nfil] [ncol] [nhilos]\n", argv[0]);
        return 1;
    }

    matriz.nfil = (unsigned int)atoi(argv[1]);
    matriz.ncol = (unsigned int)atoi(argv[2]);
    nhilos = (unsigned int)atoi(argv[3]);
  
    if(matriz.nfil < 1 || matriz.ncol < 1 || nhilos < 1){
        printf("Error de uso: nfil y ncol deben ser mayor que 0\n", argv[0]);
        return 2;
    }

    /*
        RESERVA DE MEMORIA PARA MATRICES
    */
    matriz.actual = (float **)malloc(sizeof(float *) * matriz.nfil);
    matriz.actual[0] = malloc(sizeof(float) * matriz.nfil * matriz.ncol);

    matriz.anterior = (float **)malloc(sizeof(float *) * matriz.nfil);
    matriz.anterior[0] = malloc(sizeof(float) * matriz.nfil * matriz.ncol);

    for(i = 1; i < matriz.nfil; i++){
        matriz.actual[i] = matriz.actual[i-1] + matriz.ncol;
        matriz.anterior[i] = matriz.anterior[i-1] + matriz.ncol;
    }
    
    /* 
        INICIALIZACIÓN DE TEMPERATURAS
    */
    
    for(int i = 0; i < matriz.nfil; i++){
        for(int j = 0; j < matriz.ncol; j++){
            matriz.actual[i][j] = 20;
        }
    }
    matriz.actual[0][0] = 0;
    matriz.actual[matriz.nfil-1][matriz.ncol-1] = 100;


    /*
        EJECUCIÓN 
    */
    gettimeofday(&start, NULL);
    

    hilo_func(NULL);
    gettimeofday(&end, NULL);
    status.tiempo = (end.tv_sec - start.tv_sec) +
                (end.tv_usec - start.tv_usec) / 1e6;

    /*
        IMPRIMIR RESULTADOS
    */
    #ifdef DEBUG
        for(i = 0; i < matriz.nfil; i++){
            for(id = 0; id < matriz.ncol; id++){
                printf("%f ", matriz.actual[i][id]);
            }
            printf("\n");
        }
    #endif

    /*
        LIBERACIÓN DE RECURSOS
    */
    free(matriz.actual[0]);
    free(matriz.anterior[0]);
    free(matriz.actual);
    free(matriz.anterior);

    /*
        Almacenar resultados en un fichero
    */
    f = fopen(LOG_FILE, "a");  
    if (!f) {
        perror("fopen");
        return 3;
    }
    fprintf(f, "%d,%u,%u,%f,%ld\n", nhilos, matriz.nfil, matriz.ncol, status.tiempo, status.iteraciones_temporales);

    fclose(f);
    return 0;
}


void* hilo_func(void *info) {
    float **aux;
    unsigned int i, j;
    struct Temperatura {
        float arriba, abajo, izquierda, derecha;
    };
    short int terminar = 0;

    #pragma omp parallel num_threads(nhilos) \
            private(i, j) shared(terminar, dt)
    {
        struct Temperatura t;

        while(1) {
            /*
                CÁLCULO DE TEMPERATURAS
            */
            #pragma omp single
            {
                if (dt < ESTACIONARIO) {
                    terminar = 1;
                } else {
                    aux = matriz.actual;
                    matriz.actual = matriz.anterior;
                    matriz.anterior = aux;

                    status.iteraciones_temporales++;
                }
                dt = 0.0f;
            }

            if (terminar) break;

            #pragma omp for reduction(+:dt)
            for (i = 0; i < matriz.nfil; i++){
                for(j = 0; j < matriz.ncol; j++){
                    t.arriba = i > 0 ? matriz.anterior[i][j] - matriz.anterior[i-1][j] : 0;
                    t.abajo = i < matriz.nfil - 1 ? matriz.anterior[i][j] - matriz.anterior[i+1][j] : 0;
                    t.izquierda = j > 0 ? matriz.anterior[i][j] - matriz.anterior[i][j-1] : 0;
                    t.derecha = j < matriz.ncol - 1 ? matriz.anterior[i][j] - matriz.anterior[i][j+1] : 0;
                    
                    matriz.actual[i][j] = matriz.anterior[i][j]
                    - R * t.arriba
                    - R * t.abajo
                    - R * t.izquierda
                    - R * t.derecha;  

                    if(matriz.actual[i][j] - matriz.anterior[i][j] > 0){
                        dt += matriz.actual[i][j] - matriz.anterior[i][j];
                    } else {
                        dt += matriz.anterior[i][j] - matriz.actual[i][j];
                    }   
                }
            }
        }
    }
    return NULL;
}

