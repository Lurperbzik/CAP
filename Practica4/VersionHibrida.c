#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>
#include <omp.h> 

#define R 0.1
#define ESTACIONARIO 1e-6
#define LOG_FILE "resultados.csv"

struct MatrizLocal {
    unsigned int nfil, ncol;
    float **matriz_actual;
    float **matriz_anterior;
};
struct MatrizGlobal {
    unsigned int nfil, ncol;
    float **matriz;
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
short int nhilos; // número de hilos por proceso

int main(int argc, char **argv){
    int id, nprocs, i, j;
    struct Topologia topologia;
    struct Tareas send_root;
    struct MatrizGlobal global;
    struct MatrizLocal local = {0};
    struct timeval start, end;
    float *recv_buf;
    FILE *f;
    int provided;
    
    /*
        LECTURA DE ARGUMENTOS
    */
    if (argc < 4) {
        printf("Error de uso: %s [nfil] [ncol] [nhilos]\n", argv[0]);
        return 1;
    }    
    global.nfil = (unsigned int)atoi(argv[1]);
    global.ncol = (unsigned int)atoi(argv[2]);
    nhilos = (unsigned int)atoi(argv[3]);

    if(global.nfil < 1 || global.ncol < 1){
        printf("Error de uso: nfil y ncol deben ser mayor que 0\n");
        return 2;
    }

    /*
        INICIALIZACIÓN DE MPI
    */
    gettimeofday(&start, NULL);

    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED) {
        fprintf(stderr, "MPI no soporta MPI_THREAD_FUNNELED\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    
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
        global.matriz[0] =  malloc(sizeof(float) * global.nfil * global.ncol);

        for(i = 1; i < global.nfil; i++){
            global.matriz[i] = global.matriz[i-1] + global.ncol;
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
    local.matriz_actual[0] = malloc(sizeof(float) * (local.nfil) * local.ncol);

    local.matriz_anterior = (float **)malloc(sizeof(float *) * (local.nfil));
    local.matriz_anterior[0] = malloc(sizeof(float) * (local.nfil) * local.ncol);

    for(i = 1; i < local.nfil; i++){
        local.matriz_actual[i] = local.matriz_actual[i-1] + local.ncol;
        local.matriz_anterior[i] = local.matriz_anterior[i-1] + local.ncol;
    }
 
    /* 
        ENVÍO DE TAREAS 
    */
    MPI_Type_contiguous(global.ncol, MPI_FLOAT, &send_root.tipo_fila);
    MPI_Type_commit(&send_root.tipo_fila);

    MPI_Scatterv(
        id == 0 ? global.matriz[0]  : NULL,
        id == 0 ? send_root.contadores_env  : NULL,
        id == 0 ? send_root.despls          : NULL,
        id == 0 ? send_root.tipo_fila       : MPI_DATATYPE_NULL,
        local.matriz_anterior[0],
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
        local.matriz_actual[0],
        local.nfil, 
        send_root.tipo_fila,
        id == 0 ? global.matriz[0] : NULL,
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
    free(local.matriz_actual[0]);
    free(local.matriz_anterior[0]);
    free(local.matriz_actual);
    free(local.matriz_anterior);
    free(recv_buf);
    if (id == 0){
        free(global.matriz[0]);
        free(global.matriz);
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
    struct Temperatura {
        float arriba, abajo, izquierda, derecha;
    };
    float global_dt, local_dt;
    int send_recv_count[2] = {1, 1};     // tamaño del búfer de mensajes
    MPI_Datatype send_recv_type[2] = {tipo_fila, tipo_fila}; // tipos del búfer de mensajes
    size_t send_despl[2] = {0, sizeof(float) * local.ncol * (local.nfil - 1)};
    size_t recv_despl[2] = {0, sizeof(float) * local.ncol}; 
    float *fila_superior = recv_buf;
    float *fila_inferior = recv_buf + local.ncol;

    short int terminar = 0;

    #pragma omp parallel num_threads(nhilos) private(i, j) 
    {
        struct Temperatura t;
        while(1){
            /*
                RECEPCION Y ENVIO DE/A VECINOS
            */
            #pragma omp master 
            {
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
                local_dt = 0.0f;
                status.iteraciones_temporales++;
            }
            #pragma omp barrier

            /*
                CÁLCULO DE TEMPERATURAS
            */
            #pragma omp for reduction(+:local_dt)
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
            #pragma omp master 
            {
                MPI_Allreduce(
                    &local_dt, 
                    &global_dt, 
                    1, 
                    MPI_FLOAT, 
                    MPI_SUM, 
                    MPI_COMM_WORLD
                ); 
                
                if (global_dt < ESTACIONARIO){
                    terminar = 1;
                }
                /*
                INTERCAMBIO DE PUNTEROS
                */
                float **temp = local.matriz_actual;
                local.matriz_actual = local.matriz_anterior;
                local.matriz_anterior = temp;
            }
            #pragma omp barrier 
            if (terminar) break;
        }
    }
}

/*
export OMP_NUM_THREADS=2
mpirun -np 4 ./calor 2000 2000 2

*/