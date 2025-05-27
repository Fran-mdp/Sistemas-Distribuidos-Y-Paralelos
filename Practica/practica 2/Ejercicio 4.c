#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>
#include <time.h>

int *vector;          // Vector global de enteros
int N, T;             // Tamaño del vector y número de hilos
double suma_total = 0.0; // Acumulador global de la suma
pthread_mutex_t mutex;   // Mutex para proteger la suma

// Función para obtener tiempo en segundos
double dwalltime(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec/1000000.0;
}

// Función que ejecuta cada hilo
void *sumar_parcial(void *arg){
    int id = *(int *)arg;
    int desde = (N / T) * id;
    int hasta = (id == T - 1) ? N : desde + (N / T);

    int suma_local = 0;

    // Suma local
    for(int i = desde; i < hasta; i++){
        suma_local += vector[i];
    }

    // Suma protegida con mutex (sección crítica)
    pthread_mutex_lock(&mutex);
    suma_total += suma_local;
    pthread_mutex_unlock(&mutex);

    pthread_exit(NULL);
}

// Carga el vector con números aleatorios
void inicializarVector(){
    srand(time(NULL));
    for(int i = 0; i < N; i++){
        vector[i] = rand() % 100;  // Valores entre 0 y 99
    }
}

int main(int argc, char *argv[]){
    if(argc != 3){
        printf("Uso: %s <nro_hilos> <dimension_N>\n", argv[0]);
        return -1;
    }

    T = atoi(argv[1]);  // Número de hilos
    N = atoi(argv[2]);  // Tamaño del vector

    vector = malloc(sizeof(int) * N);
    pthread_t threads[T];
    int ids[T];
    pthread_mutex_init(&mutex, NULL);

    inicializarVector();  // Llenar el vector con valores aleatorios

    double inicio = dwalltime();

    // Crear hilos
    for(int i = 0; i < T; i++){
        ids[i] = i;
        pthread_create(&threads[i], NULL, sumar_parcial, (void *)&ids[i]);
    }

    // Esperar que todos los hilos terminen
    for(int i = 0; i < T; i++){
        pthread_join(threads[i], NULL);
    }

    double fin = dwalltime();

    // Calcular promedio
    double promedio = suma_total / N;

    // Mostrar resultados
    printf("Promedio: %.2f\n", promedio);
    printf("Tiempo con %d hilos y %d elementos: %.5f segundos\n", T, N, fin - inicio);

    free(vector);
    pthread_mutex_destroy(&mutex);
    return 0;
}
