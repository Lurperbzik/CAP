#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define R 0.1
#define ESTACIONARIO 1e-6
#define LOG_FILE "resultados.csv"

void* hilo_func(void *);

struct Matriz {
    unsigned int nfil, ncol;
    float **actual;
    float **anterior;
};
struct Matriz matriz = {0};

float dt = 0;

struct Status {
    double tiempo;
    long iteraciones_temporales;
};
struct Status status = {0};

int main(int argc, char * argv[]) {
    /*
        VARIABLES LOCALES
    */
    float *datos_actual, *datos_anterior;
    int i;
    long id;
    struct timeval start, end;
    FILE *f;
    

    /*
        LECTURA DE ARGUMENTOS
    */
    if (argc < 3) {
        printf("Error de uso: %s [nfil] [ncol]\n", argv[0]);
        return 1;
    }

    matriz.nfil = (unsigned int)atoi(argv[1]);
    matriz.ncol = (unsigned int)atoi(argv[2]);
  
    if(matriz.nfil < 1 || matriz.ncol < 1){
        printf("Error de uso: nfil y ncol deben ser mayor que 0\n", argv[0]);
        return 2;
    }

    /*
        RESERVA DE MEMORIA PARA MATRICES
    */
    matriz.actual = (float **)malloc(sizeof(float *) * matriz.nfil);
    datos_actual = malloc(sizeof(float) * matriz.nfil * matriz.ncol);

    matriz.anterior = (float **)malloc(sizeof(float *) * matriz.nfil);
    datos_anterior = malloc(sizeof(float) * matriz.nfil * matriz.ncol);
    
    for(i = 0; i < matriz.nfil; i++){
        matriz.actual[i] = datos_actual + i * matriz.ncol;
        matriz.anterior[i] = datos_anterior + i * matriz.ncol;
    }
    
    /* 
        INICIALIZACIÓN DE TEMPERATURAS
    */
    
    for(int i = 0; i < matriz.nfil; i++){
        for(int j = 0; j < matriz.ncol; j++){
            matriz.anterior[i][j] = 20;
        }
    }
    matriz.anterior[0][0] = 0;
    matriz.anterior[matriz.nfil-1][matriz.ncol-1] = 100;

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
    free(datos_actual);
    free(datos_anterior);
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
    fprintf(f, "%d,%u,%u,%f,%ld\n", 1, matriz.nfil, matriz.ncol, status.tiempo, status.iteraciones_temporales);

    fclose(f);
    return 0;
}


void* hilo_func(void *info) {
    float **aux;
    unsigned int i, j;
    struct Temperatura {
        float arriba, abajo, izquierda, derecha;
    };
    struct Temperatura t;

    while(1) {
        /*
            CÁLCULO DE TEMPERATURAS
        */
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
        
        status.iteraciones_temporales++;
    
        // intercambiar punteros de matrices 
        aux = matriz.actual;
        matriz.actual = matriz.anterior;
        matriz.anterior = aux;

        if (dt < ESTACIONARIO) {
            return NULL;
        }
        
        dt = 0;
    }
}

