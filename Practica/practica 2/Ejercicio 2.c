#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>
#include <time.h>

int *vector;            // Vector global
int N, T, X;            // Tamaño del vector, cantidad de hilos, valor a buscar
int ocurrencias = 0;    // Resultado final
pthread_mutex_t mutex;  // Para sincronizar el acceso a la variable global

// Función para medir tiempo real
double dwalltime(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec/1000000.0;
}

// Función que ejecuta cada hilo
void *contar(void *arg){
    int id = *(int *)arg;
    int desde = (N / T) * id;
    int hasta = (id == T - 1) ? N : desde + (N / T);
    int cuenta_local = 0;

    // Cada hilo cuenta sus ocurrencias localmente
    for(int i = desde; i < hasta; i++){
        if(vector[i] == X){
            cuenta_local++;
        }
    }

    // Suma los resultados al contador global (con mutex para evitar condiciones de carrera)
    pthread_mutex_lock(&mutex);
    ocurrencias += cuenta_local;
    pthread_mutex_unlock(&mutex);

    pthread_exit(NULL);
}

// Inicializa el vector con valores aleatorios
void inicializarVector(){
    srand(time(NULL));
    for(int i = 0; i < N; i++){
        vector[i] = rand() % 10;  // Números del 0 al 9
    }
}

int main(int argc, char *argv[]){
    if(argc != 4){
        printf("Uso: %s <nro_hilos> <dimension_N> <valor_X>\n", argv[0]);
        return -1;
    }

    T = atoi(argv[1]);
    N = atoi(argv[2]);
    X = atoi(argv[3]);

    pthread_t threads[T];
    int ids[T];
    vector = malloc(sizeof(int) * N);
    pthread_mutex_init(&mutex, NULL);

    inicializarVector();

    double inicio = dwalltime();

    // Crear hilos
    for(int i = 0; i < T; i++){
        ids[i] = i;
        pthread_create(&threads[i], NULL, contar, (void *)&ids[i]);
    }

    // Esperar a que todos terminen
    for(int i = 0; i < T; i++){
        pthread_join(threads[i], NULL);
    }

    double fin = dwalltime();

    // Resultado final
    printf("Cantidad de ocurrencias de %d: %d\n", X, ocurrencias);
    printf("Tiempo con %d hilos y vector de %d elementos: %.5f segundos\n", T, N, fin - inicio);

    free(vector);
    pthread_mutex_destroy(&mutex);

    return 0;
}
