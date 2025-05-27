#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>
#include <limits.h>
#include <time.h>

int *vector;
int N, T;

int min_global, max_global;
pthread_mutex_t mutex;

// Estructura para guardar min y max local de cada hilo
typedef struct {
    int min;
    int max;
} Resultado;

double dwalltime(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec/1000000.0;
}

// Función que ejecuta cada hilo
void *buscar_min_max(void *arg){
    int id = *(int *)arg;
    int desde = (N / T) * id;
    int hasta = (id == T - 1) ? N : desde + (N / T);

    int min_local = INT_MAX;
    int max_local = INT_MIN;

    // Buscar min y max en la porción del hilo
    for(int i = desde; i < hasta; i++){
        if(vector[i] < min_local) min_local = vector[i];
        if(vector[i] > max_local) max_local = vector[i];
    }

    // Proteger actualizaciones globales con mutex
    pthread_mutex_lock(&mutex);
    if(min_local < min_global) min_global = min_local;
    if(max_local > max_global) max_global = max_local;
    pthread_mutex_unlock(&mutex);

    pthread_exit(NULL);
}

// Genera números aleatorios en el vector
void inicializarVector(){
    srand(time(NULL));
    for(int i = 0; i < N; i++){
        vector[i] = rand() % 10000;  // valores entre 0 y 9999
    }
}

int main(int argc, char *argv[]){
    if(argc != 3){
        printf("Uso: %s <nro_hilos> <dimension_N>\n", argv[0]);
        return -1;
    }

    T = atoi(argv[1]);
    N = atoi(argv[2]);

    vector = malloc(sizeof(int) * N);
    pthread_t threads[T];
    int ids[T];
    pthread_mutex_init(&mutex, NULL);

    inicializarVector();

    min_global = INT_MAX;
    max_global = INT_MIN;

    double inicio = dwalltime();

    for(int i = 0; i < T; i++){
        ids[i] = i;
        pthread_create(&threads[i], NULL, buscar_min_max, (void *)&ids[i]);
    }

    for(int i = 0; i < T; i++){
        pthread_join(threads[i], NULL);
    }

    double fin = dwalltime();

    printf("Mínimo global: %d\n", min_global);
    printf("Máximo global: %d\n", max_global);
    printf("Tiempo con %d hilos y vector de %d elementos: %.5f segundos\n", T, N, fin - inicio);

    free(vector);
    pthread_mutex_destroy(&mutex);
    return 0;
}
