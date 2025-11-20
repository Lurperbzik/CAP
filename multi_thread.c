#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
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

struct Barrier {
    pthread_mutex_t mutex; 
    pthread_cond_t cond;
    unsigned short int cont;
};
struct Barrier barrier = {0};

sem_t sem_dt;
unsigned short int nhilos;
float dt = 0;
char fin_programa = -1;

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
    pthread_t *thread_handles;
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
    nhilos = atoi(argv[3]);

    if(matriz.nfil < 1 || matriz.ncol < 1 || nhilos < 1){
        printf("Error de uso: nfil, ncol y nhilos deben ser mayor que 0\n", argv[0]);
        return 2;
    }

    /*
        INICIALIZACIÓN DE RECURSOS DE SINCRONIZACIÓN
    */
    sem_init(&sem_dt, 0, 1);
    pthread_mutex_init(&barrier.mutex, NULL);
    pthread_cond_init(&barrier.cond, NULL);

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
        CREACIÓN DE HILOS
    */

    thread_handles = malloc(nhilos*sizeof(pthread_t));

    gettimeofday(&start, NULL);

    for (id = 0; id < nhilos; id++){
        if(pthread_create(&thread_handles[id], NULL, hilo_func, (void*)id)){
            printf("Error creando hilo\n");
        } 
    }

    for (i = 0; i < nhilos; i++){
       pthread_join(thread_handles[i], NULL);
    }
   
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
    free(thread_handles);
   
    sem_destroy(&sem_dt);
    pthread_mutex_destroy(&barrier.mutex);
    pthread_cond_destroy(&barrier.cond);
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
    fprintf(f, "%hu,%u,%u,%f,%ld\n", nhilos, matriz.nfil, matriz.ncol, status.tiempo, status.iteraciones_temporales);

    fclose(f);
    return 0;
   
}


void* hilo_func(void *info) {
    long id = (long)info;
    float **aux;
    /*
        CÁLCULO DE FILAS
    */
    unsigned int nfilas = matriz.nfil / nhilos;
    unsigned int primera_fila, ultima_fila;
    unsigned int i, j; 
    float local_dt; 
    struct Temperatura {
        float arriba, abajo, izquierda, derecha;
    };
    struct Temperatura t;

    /*
        ASIGNACIÓN DE TAREAS A CADA HILO
    */

    if (matriz.nfil % nhilos != 0){
        if (id < matriz.nfil%nhilos){
            primera_fila = id * nfilas + id;
            ultima_fila = primera_fila + nfilas;
        } else {
            primera_fila = id * nfilas + matriz.nfil%nhilos;
            ultima_fila = primera_fila + nfilas - 1;
        }
    } else {
        primera_fila = id * nfilas;
        ultima_fila = primera_fila + nfilas - 1;
    }

    while(1) {
        /*
            CÁLCULO DE TEMPERATURAS
        */
        local_dt = 0;
        for (i = primera_fila; i <= ultima_fila; i++){
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
                    local_dt += matriz.actual[i][j] - matriz.anterior[i][j];
                } else {
                    local_dt += matriz.anterior[i][j] - matriz.actual[i][j];
                }   
            }
        } 
        sem_wait(&sem_dt);
        dt += local_dt;
        sem_post(&sem_dt);
        /*
            BARRERA DE SINCRONIZACIÓN
        */
        pthread_mutex_lock(&barrier.mutex);
        barrier.cont++;
        if(barrier.cont == nhilos){
            // intercambiar punteros de matrices 
            aux = matriz.actual;
            matriz.actual = matriz.anterior;
            matriz.anterior = aux;
            
            status.iteraciones_temporales++;

            /*
                COMPRUEBA SI HAY ESTADO DE TEMPERATURA ESTACIONARIO
            */
            if (dt < ESTACIONARIO) {
                fin_programa = 0;
                pthread_cond_broadcast(&barrier.cond);
                pthread_mutex_unlock(&barrier.mutex); 
                return NULL;
            }
            barrier.cont = 0;
            dt = 0;
            pthread_cond_broadcast(&barrier.cond);
        } else {
            // Desbloque la barrara y bloquea la ejecución de el hilo
            // Se pone el wait en un while por si otro evento despierta el hilo
            while(pthread_cond_wait(&barrier.cond, &barrier.mutex) != 0);
            /*
                COMPRUEBA SI HAY ESTADO DE TEMPERATURA ESTACIONARIO
            */
            if(fin_programa == 0){
                pthread_mutex_unlock(&barrier.mutex); 
                return NULL;
            }
        }
        pthread_mutex_unlock(&barrier.mutex);
    }
}


