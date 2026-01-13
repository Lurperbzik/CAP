#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>

#define R 0.1
#define ESTACIONARIO 1e-6
#define LOG_FILE "resultados.csv"

struct MatrizLocal {
    unsigned int nfil, ncol;
    float **matriz_actual;
    float **matriz_anterior;
    float *datos_actual;
    float *datos_anterior; 
};
struct MatrizGlobal {
    unsigned int nfil, ncol;
    float **matriz;
    float *datos;
};
struct Topologia {
    int dims[1];
    int periods[1];
};
struct Tareas {
    int *contadores_env;
    int *despls;
    MPI_Datatype tipo_fila;
};
struct Status {
    double tiempo;
    long iteraciones_temporales;
};

void computar(struct MatrizLocal, MPI_Comm, float *, int, int, MPI_Datatype);
struct Status status = {0};


int main(int argc, char **argv){
    int id, nprocs, i, j;
    struct Topologia topologia;
    struct Tareas send_root;
    struct MatrizGlobal global;
    struct MatrizLocal local = {0};
    struct timeval start, end;
    float *recv_buf;
    FILE *f;
    
    /*
        LECTURA DE ARGUMENTOS
    */
    if (argc < 3) {
        printf("Error de uso: %s [nfil] [ncol]\n", argv[0]);
        return 1;
    }    
    global.nfil = (unsigned int)atoi(argv[1]);
    global.ncol = (unsigned int)atoi(argv[2]);

    if(global.nfil < 1 || global.ncol < 1){
        printf("Error de uso: nfil y ncol deben ser mayor que 0\n");
        return 2;
    }

    /*
        INICIALIZACIÓN DE MPI
    */
    gettimeofday(&start, NULL);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    /*
        PROCESO 0
    */
    if(id == 0){
        /*
            RESERVA DE MEMORIA 
        */
        global.matriz = (float **)malloc(sizeof(float *) * (global.nfil));
        global.datos = malloc(sizeof(float) * (global.nfil) * global.ncol);
        for(i = 0; i < global.nfil; i++){
            global.matriz[i] = global.datos + i * global.ncol;
        }

        /* 
            INICIALIZACIÓN DE TEMPERATURAS
        */
        for(i = 0; i < global.nfil; i++){
            for(j = 0; j < global.ncol; j++){
               global.matriz[i][j] = 20;
            }
        }
        global.matriz[0][0] = 0;
        global.matriz[global.nfil-1][global.ncol-1] = 100;   

        /*
            DATOS DE SCATTERV
        */
        send_root.contadores_env = malloc(sizeof(int) * nprocs);
        send_root.despls = malloc(sizeof(int) * nprocs);

        for(i = 0; i < nprocs; i++){
            send_root.contadores_env[i] = (i < global.nfil % nprocs)
                ? global.nfil / nprocs + 1 
                : global.nfil / nprocs;
            
            send_root.despls[i] = (i == 0) 
                ? 0 :
                send_root.despls[i-1] + send_root.contadores_env[i-1];
        }
    }
      

    /*
        ASIGNACIÓN DE TAREAS A CADA PROCESO
    */
    local.nfil = (id < global.nfil % nprocs)
        ? global.nfil / nprocs + 1 
        : global.nfil / nprocs;
    local.ncol = global.ncol;

    /*
        RESERVA DE MEMORIA PARA MATRICES LOCALES
    */
    local.matriz_actual = (float **)malloc(sizeof(float *) * (local.nfil));
    local.datos_actual = malloc(sizeof(float) * (local.nfil) * local.ncol);

    local.matriz_anterior = (float **)malloc(sizeof(float *) * (local.nfil));
    local.datos_anterior = malloc(sizeof(float) * (local.nfil) * local.ncol);

    for(i = 0; i < local.nfil; i++){
        local.matriz_actual[i] = local.datos_actual + i * local.ncol;
        local.matriz_anterior[i] = local.datos_anterior + i * local.ncol;
    }
    /*
        REGISTRAR TIEMPO DE INICIO
    */
 
    /* 
        ENVÍO DE TAREAS 
    */
    MPI_Type_contiguous(global.ncol, MPI_FLOAT, &send_root.tipo_fila);
    MPI_Type_commit(&send_root.tipo_fila);

    MPI_Scatterv(
        id == 0 ? global.datos  : NULL,
        id == 0 ? send_root.contadores_env  : NULL,
        id == 0 ? send_root.despls          : NULL,
        id == 0 ? send_root.tipo_fila       : MPI_DATATYPE_NULL,
        local.datos_anterior,
        local.nfil,
        send_root.tipo_fila,
        0,
        MPI_COMM_WORLD
    );

    /* RESERVA DE MEMORIA PARA BÚFERES DE RECEPCIÓN */
    recv_buf = malloc(local.ncol * 2 * sizeof(float));
    
    /* 
        CREACIÓN TOPOLOGÍA CARTESIANA 1D -> TOPOLOGÍA DE VECINOS
    */
    topologia.dims[0] = nprocs;
    topologia.periods[0] = 0;
    MPI_Comm comm1D;
    MPI_Cart_create(
        MPI_COMM_WORLD, 
        1, 
        topologia.dims, 
        topologia.periods, 
        0, 
        &comm1D
    );

    /*
        CÁLCULOS
    */
    computar(local, comm1D, recv_buf, id, nprocs, send_root.tipo_fila);

    /*
        RECEPCIÓN DE DATOS EN PROCESO RAÍZ
    */

    MPI_Gatherv(
        local.datos_actual,
        local.nfil, 
        send_root.tipo_fila,
        id == 0 ? global.datos : NULL,
        id == 0 ? send_root.contadores_env : NULL,
        id == 0 ? send_root.despls : NULL,
        id == 0 ? send_root.tipo_fila : MPI_DATATYPE_NULL,
        0,
        MPI_COMM_WORLD
    );
    /*
        REGISTRAR TIEMPO DE FIN
    */
    if(id == 0){
        gettimeofday(&end, NULL);
        status.tiempo = (end.tv_sec - start.tv_sec) +
            (end.tv_usec - start.tv_usec) / 1e6;
    }
    
    #ifdef DEBUG
        if(id == 0){
            for(i = 0; i < global.nfil; i++){
                for(j = 0; j < global.ncol; j++){
                    printf("%f ", global.matriz[i][j]);
                }
                printf("\n");
            }
            printf("Iteraciones : %d\n", status.iteraciones_temporales);
        }
    #endif

    /*
        LIBERAR RECURSOS
    */  
    free(local.datos_actual);
    free(local.datos_anterior);
    free(local.matriz_actual);
    free(local.matriz_anterior);
    free(recv_buf);
    if (id == 0){
        free(global.matriz);
        free(global.datos);
    }
    /*
        ALMACENAR LOS DATOS EN UN FICHERO
    */
    if(id == 0){
        f = fopen(LOG_FILE, "a");  
        if (!f) {
            perror("fopen");
            return 3;
        }
        fprintf(f, "%hu,%u,%u,%f,%ld\n", nprocs, global.nfil, global.ncol, status.tiempo, status.iteraciones_temporales);
        fclose(f);
   }
    /*
        FIN
    */
    MPI_Finalize();
    return 0;
}


void computar(struct MatrizLocal local, MPI_Comm comm1D, float *recv_buf, int id, int nprocs, MPI_Datatype tipo_fila){
    int i, j;
    float **temp;
    struct Temperatura {
        float arriba, abajo, izquierda, derecha;
    };
    struct Temperatura t;
    float local_dt, global_dt;

    int send_recv_count[2] = {1, 1};     // tamaño del búfer de mensajes
    MPI_Datatype send_recv_type[2] = {tipo_fila, tipo_fila}; // tipos del búfer de mensajes
    size_t send_despl[2] = {0, sizeof(float) * local.ncol * (local.nfil - 1)};
    size_t recv_despl[2] = {0, sizeof(float) * local.ncol}; 
    float *fila_superior = recv_buf;
    float *fila_inferior = recv_buf + local.ncol;
    
    while(1){
        local_dt = 0;
        status.iteraciones_temporales++;
        /*
            RECEPCION Y ENVIO DE/A VECINOS
        */
        MPI_Neighbor_alltoallw(
            local.matriz_anterior[0], 
            send_recv_count, 
            send_despl, 
            send_recv_type,
            recv_buf, 
            send_recv_count, 
            recv_despl, 
            send_recv_type,
            comm1D
        );

        /*
            CÁLCULO DE TEMPERATURAS
        */
        for(i = 0; i < local.nfil; i++){
            for(j = 0; j < local.ncol; j++){
                
                if (i == 0){
                    t.arriba = id != 0 ? local.matriz_anterior[i][j] - fila_superior[j] : 0;
                } else{
                    t.arriba = local.matriz_anterior[i][j] - local.matriz_anterior[i-1][j];
                } 
                if (i == local.nfil - 1){
                    t.abajo = id != nprocs - 1 ? local.matriz_anterior[i][j] - fila_inferior[j] : 0;
                } else {
                    t.abajo = local.matriz_anterior[i][j] - local.matriz_anterior[i+1][j];
                } 

                t.izquierda = j > 0 ? local.matriz_anterior[i][j] - local.matriz_anterior[i][j-1] : 0;
                t.derecha = j < local.ncol - 1 ? local.matriz_anterior[i][j] - local.matriz_anterior[i][j+1] : 0;                
                
                local.matriz_actual[i][j] = local.matriz_anterior[i][j]
                    - R * t.arriba
                    - R * t.abajo
                    - R * t.izquierda
                    - R * t.derecha;  

                if(local.matriz_actual[i][j] - local.matriz_anterior[i][j] > 0){
                    local_dt += local.matriz_actual[i][j] - local.matriz_anterior[i][j];
                } else {
                    local_dt += local.matriz_anterior[i][j] - local.matriz_actual[i][j];
                }   
            }
        }

        /*
            ENVÍO DE DIFERENCIAL DE TEMPERATURA
        */
        MPI_Allreduce(
            &local_dt, 
            &global_dt, 
            1, 
            MPI_FLOAT, 
            MPI_SUM, 
            MPI_COMM_WORLD
        ); // también funciona como barrera de sincronización
        
        if (global_dt < ESTACIONARIO){
            return;
        } 
        /*
            INTERCAMBIO DE PUNTEROS
        */
        temp = local.matriz_actual;
        local.matriz_actual = local.matriz_anterior;
        local.matriz_anterior = temp;
    }

}