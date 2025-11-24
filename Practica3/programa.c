#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define R 0.1
#define ESTACIONARIO 1e-6
#define LOG_FILE "resultados.csv"
struct Matriz {
    unsigned int nfil, ncol;
    float **actual;
    float **anterior;
};
int iteraciones = 0;
void computar(struct Matriz, MPI_Comm, float *, int, int);

int main(int argc, char **argv){
    int id, nprocs;
    int dims[1];                // tamaño de la malla 1D
    int periods[1] = {0};       // sin periodicidad
    struct Matriz matriz = {0};
    float *datos_actual, *datos_anterior;
    float *recv_buf;
    int i, j;
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
        printf("Error de uso: nfil y ncol deben ser mayor que 0\n");
        return 2;
    }

    /*
        INICIALIZACIÓN DE MPI
    */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    /*
        ASIGNACIÓN DE TAREAS A CADA PROCESO
    */
    matriz.nfil = id < matriz.nfil % nprocs ? matriz.nfil / nprocs + 1 : matriz.nfil / nprocs;

    /*
        RESERVA DE MEMORIA PARA MATRICES
    */
    matriz.actual = (float **)malloc(sizeof(float *) * (matriz.nfil));
    datos_actual = malloc(sizeof(float) * (matriz.nfil) * matriz.ncol);

    matriz.anterior = (float **)malloc(sizeof(float *) * (matriz.nfil));
    datos_anterior = malloc(sizeof(float) * (matriz.nfil) * matriz.ncol);
    
    for(i = 0; i < matriz.nfil; i++){
        matriz.actual[i] = datos_actual + i * matriz.ncol;
        matriz.anterior[i] = datos_anterior + i * matriz.ncol;
    }

    recv_buf = malloc(matriz.ncol * 2 * sizeof(float));
    
    /* 
        INICIALIZACIÓN DE TEMPERATURAS
    */
    for(i = 0; i < matriz.nfil; i++){
        for(j = 0; j < matriz.ncol; j++){
            matriz.anterior[i][j] = 20;
        }
    }
    if (id == 0) {
        matriz.anterior[0][0] = 0;
    }
    if(id == nprocs - 1){
        matriz.anterior[matriz.nfil-1][matriz.ncol-1] = 100;
    }

    /* 
        CREACIÓN TOPOLOGÍA CARTESIANA 1D -> TOPOLOGÍA DE VECINOS
    */
    dims[0] = nprocs;
    MPI_Comm comm1D;
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 0, &comm1D);

    /*
        CÁLCULOS
    */
    computar(matriz, comm1D, recv_buf, id, nprocs);

    for(i = 0; i < matriz.nfil; i++){
        for(j = 0; j < matriz.ncol; j++){
            printf("%f ", matriz.actual[i][j]);
        }
        printf("\n");
    }
    printf("Iteraciones : %d\n", iteraciones);

    /*
        FIN
    */
    MPI_Finalize();
    return 0;
}

void computar(struct Matriz matriz, MPI_Comm comm1D, float *recv_buf, int id, int nprocs){
    int i, j;
    float **temp;
    struct Temperatura {
        float arriba, abajo, izquierda, derecha;
    };
    struct Temperatura t;
    float local_dt, global_dt;

    int send_recv_count[2] = {matriz.ncol, matriz.ncol};     // tamaño del búfer de mensajes
    MPI_Datatype send_recv_type[2] = {MPI_FLOAT, MPI_FLOAT}; // tipos del búfer de mensajes
     // desplazamientos locales
    size_t send_despl[2] = {0, sizeof(float) * matriz.ncol * (matriz.nfil - 1)};
    size_t recv_despl[2] = {0, sizeof(float) * matriz.ncol}; 
    float *fila_superior = recv_buf;
    float *fila_inferior = recv_buf + matriz.ncol;
    
    while(1){
        local_dt = 0;
        iteraciones++;
        /*
            RECEPCION Y ENVIO DE/A VECINOS
        */
        MPI_Neighbor_alltoallw(
            matriz.anterior[0], send_recv_count, send_despl, send_recv_type,
            recv_buf, send_recv_count, recv_despl, send_recv_type,
            comm1D
        );

        /*
            CÁLCULO DE TEMPERATURAS
        */
        for(i = 0; i < matriz.nfil; i++){
            for(j = 0; j < matriz.ncol; j++){
                
                if (i == 0){
                    t.arriba = id != 0 ? matriz.anterior[i][j] - fila_superior[j] : 0;
                } else{
                    t.arriba = matriz.anterior[i][j] - matriz.anterior[i-1][j];
                } 
                if (i == matriz.nfil - 1){
                    t.abajo = id != nprocs - 1 ? matriz.anterior[i][j] - fila_inferior[j] : 0;
                } else {
                    t.abajo = matriz.anterior[i][j] - matriz.anterior[i+1][j];
                } 

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

        /*
            ENVÍO DE DIFERENCIAL DE TEMPERATURA
        */
        MPI_Allreduce(&local_dt, &global_dt, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); // también funciona como barrera de sincronización
        if (global_dt < ESTACIONARIO){
            return;
        } 

        /*
            INTERCAMBIO DE PUNTEROS
        */
        temp = matriz.actual;
        matriz.actual = matriz.anterior;
        matriz.anterior = temp;
    }

}