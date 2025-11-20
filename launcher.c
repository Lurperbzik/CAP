#include <stdio.h>
#include <spawn.h>
#include <sys/wait.h>

extern char **environ;

int main(void) {

    char *nhilos[] = { "2", "4", "6", "8", "10", "12", "14", "16", "18", "20" };
    char *matriz_dim[] = { "100", "500", "750", "1000", "2000" };
    int status;
    int i, j, k;
    pid_t pid;
    char *argv[] = {
        "./single_thread",
        matriz_dim[0],
        matriz_dim[0],
        NULL
    };

    char *argv2[] = {
        "./multi_thread",
        matriz_dim[0],
        matriz_dim[0],
        nhilos[0],
        NULL
    };

    for(i = 0; i < 5; i++){
        argv[1] = matriz_dim[i];
        argv[2] = matriz_dim[i]; 
        for(j = 0; j < 5; j++){
            if (posix_spawn(&pid, "./single_thread", NULL, NULL, argv, environ) != 0) {
                perror("posix_spawn");
            }
            waitpid(pid, &status, 0);
            if(status != 0){
                perror("single_thread");
            }
        
        }
    }

    for (k = 0; k < 10; k++){
        argv2[3] = nhilos[k];
        for(i = 0; i < 5; i++){
            argv2[1] = matriz_dim[i];
            argv2[2] = matriz_dim[i];
            for(j = 0; j < 5; j++){
                if (posix_spawn(&pid, "./multi_thread", NULL, NULL, argv2, environ) != 0) {
                    perror("posix_spawn");
                }
                waitpid(pid, &status, 0);
                if(status != 0){
                    perror("single_thread");
                }
            
            }
        }
    }
 
    
    return 0;
}
