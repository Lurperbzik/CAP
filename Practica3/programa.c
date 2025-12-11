#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>

// MPI_Gatherv para recibir
// MPI_Scatterv para enviar

#define R 0.1
#define ESTACIONARIO 1e-6
#define LOG_FILE "resultados.csv"
struct Matriz {
    unsigned int nfil, ncol;
    float **actual;
    float **anterior;
};
void computar(struct Matriz, MPI_Comm, float *, int, int, MPI_Datatype);
struct Status {
    double tiempo;
    long iteraciones_temporales;
};
struct Status status = {0};

int main(int argc, char **argv){
    int id, nprocs;
    int dims[1];                // tamaño de la malla 1D
    int periods[1] = {0};       // sin periodicidad
    struct Matriz matriz = {0};
    float **matriz_completa;
    unsigned int fil_completas, col_completas;
    float *datos_actual, *datos_anterior, *datos_completa;
    float *recv_buf;
    int i, j;
    struct timeval start, end;
    FILE *f;
    // DATOS PARA SCATTERV / GATHERV
    int *contadores_env; // vector con el número de elementos por proceso
    int *despls; // Desplazamientos
    int contador_desp = 0;
    MPI_Datatype rowtype;
    /*
        LECTURA DE ARGUMENTOS
    */
    if (argc < 3) {
        printf("Error de uso: %s [nfil] [ncol]\n", argv[0]);
        return 1;
    }    
    fil_completas = (unsigned int)atoi(argv[1]);
    col_completas = (unsigned int)atoi(argv[2]);

    if(fil_completas < 1 || col_completas < 1){
        printf("Error de uso: nfil y ncol deben ser mayor que 0\n");
        return 2;
    }
    matriz.ncol = col_completas;

    /*
        INICIALIZACIÓN DE MPI
    */
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
        matriz_completa = (float **)malloc(sizeof(float *) * (fil_completas));
        datos_completa = malloc(sizeof(float) * (fil_completas) * col_completas);
        for(i = 0; i < fil_completas; i++){
            matriz_completa[i] = datos_completa + i * col_completas;
        }

        /* 
            INICIALIZACIÓN DE TEMPERATURAS
        */
        for(i = 0; i < fil_completas; i++){
            for(j = 0; j < col_completas; j++){
               matriz_completa[i][j] = 20;
            }
        }
        
        matriz_completa[0][0] = 0;
        matriz_completa[fil_completas-1][col_completas-1] = 100;   

        /*
            DATOS DE SCATTERV
        */
        contadores_env = malloc(sizeof(int) * nprocs);
        despls = malloc(sizeof(int) * nprocs);

        for(i = 0; i < nprocs; i++){
            contadores_env[i] = (i < fil_completas % nprocs)
                ? fil_completas / nprocs + 1 
                : fil_completas / nprocs;
            
            despls[i] = (i == 0) 
                ? 0 :
                despls[i-1] + contadores_env[i-1];
        }
    }
      

    /*
        ASIGNACIÓN DE TAREAS A CADA PROCESO
    */
    matriz.nfil = (id < fil_completas % nprocs)
        ? fil_completas / nprocs + 1 
        : fil_completas / nprocs;

    /*
        RESERVA DE MEMORIA PARA MATRICES LOCALES
    */
    matriz.actual = (float **)malloc(sizeof(float *) * (matriz.nfil));
    datos_actual = malloc(sizeof(float) * (matriz.nfil) * matriz.ncol);

    matriz.anterior = (float **)malloc(sizeof(float *) * (matriz.nfil));
    datos_anterior = malloc(sizeof(float) * (matriz.nfil) * matriz.ncol);

    for(i = 0; i < matriz.nfil; i++){
        matriz.actual[i] = datos_actual + i * matriz.ncol;
        matriz.anterior[i] = datos_anterior + i * matriz.ncol;
    }
    /*
        REGISTRAR TIEMPO DE INICIO
    */
    if(id == 0){
        gettimeofday(&start, NULL);
    }
    /* 
        ENVÍO DE TAREAS 
    */
    MPI_Type_contiguous(col_completas, MPI_FLOAT, &rowtype);
    MPI_Type_commit(&rowtype);

    MPI_Scatterv(
        id == 0 ? datos_completa  : NULL,
        id == 0 ? contadores_env  : NULL,
        id == 0 ? despls          : NULL,
        id == 0 ? rowtype         : MPI_DATATYPE_NULL,
        datos_anterior,
        matriz.nfil,
        rowtype,
        0,
        MPI_COMM_WORLD
    );

    /* RESERVA DE MEMORIA PARA BÚFERES DE RECEPCIÓN */
    recv_buf = malloc(matriz.ncol * 2 * sizeof(float));
    
    /* 
        CREACIÓN TOPOLOGÍA CARTESIANA 1D -> TOPOLOGÍA DE VECINOS
    */
    dims[0] = nprocs;
    MPI_Comm comm1D;
    MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 0, &comm1D);

    /*
        CÁLCULOS
    */
    computar(matriz, comm1D, recv_buf, id, nprocs, rowtype);

    /*
        RECEPCIÓN DE DATOS EN PROCESO RAÍZ
    */

    MPI_Gatherv(
        datos_actual,
        matriz.nfil, 
        rowtype,
        id == 0 ? datos_completa : NULL,
        id == 0 ? contadores_env : NULL,
        id == 0 ? despls : NULL,
        id == 0 ? rowtype : MPI_DATATYPE_NULL,
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
            for(i = 0; i < fil_completas; i++){
                for(j = 0; j < col_completas; j++){
                    printf("%f ", matriz_completa[i][j]);
                }
                printf("\n");
            }
            printf("Iteraciones : %d\n", status.iteraciones_temporales);
        }
    #endif

    /*
        LIBERAR RECURSOS
    */  
    free(datos_actual);
    free(datos_anterior);
    free(matriz.actual);
    free(matriz.anterior);
    free(recv_buf);
    if (id == 0){
        free(matriz_completa);
        free(datos_completa);
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
        fprintf(f, "%hu,%u,%u,%f,%ld\n", nprocs, fil_completas, col_completas, status.tiempo, status.iteraciones_temporales);
        fclose(f);
   }
    /*
        FIN
    */
    MPI_Finalize();
    return 0;
}


void computar(struct Matriz matriz, MPI_Comm comm1D, float *recv_buf, int id, int nprocs, MPI_Datatype rowtype){
    int i, j;
    float **temp;
    struct Temperatura {
        float arriba, abajo, izquierda, derecha;
    };
    struct Temperatura t;
    float local_dt, global_dt;

    int send_recv_count[2] = {1, 1};     // tamaño del búfer de mensajes
    MPI_Datatype send_recv_type[2] = {rowtype, rowtype}; // tipos del búfer de mensajes
     // desplazamientos locales
    size_t send_despl[2] = {0, sizeof(float) * matriz.ncol * (matriz.nfil - 1)};
    size_t recv_despl[2] = {0, sizeof(float) * matriz.ncol}; 
    float *fila_superior = recv_buf;
    float *fila_inferior = recv_buf + matriz.ncol;
    
    while(1){
        local_dt = 0;
        status.iteraciones_temporales++;
        /*
            RECEPCION Y ENVIO DE/A VECINOS
        */
        MPI_Neighbor_alltoallw(
            matriz.anterior[0], 
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
        temp = matriz.actual;
        matriz.actual = matriz.anterior;
        matriz.anterior = temp;
    }

}