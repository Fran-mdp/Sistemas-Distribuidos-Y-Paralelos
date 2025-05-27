#include <stdio.h>      // Para printf
#include <stdlib.h>     // Para malloc, free, atoi
#include <pthread.h>    // Para usar hilos (Pthreads)
#include <sys/time.h>   // Para medir tiempo de ejecución

// Variables globales
double **A, **B, **C;  // Matrices de entrada A, B y matriz resultado C
int N, T;              // N: tamaño de la matriz, T: número de hilos

// Función para medir el tiempo de ejecución (tiempo real)
double dwalltime(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec/1000000.0;
}

// Función que ejecuta cada hilo (recibe su ID)
void *multiplicar(void *arg){
    int id = *(int *)arg;

    // Calcula desde qué fila hasta cuál fila le toca a este hilo
    int desde = (N / T) * id;
    int hasta = (id == T - 1) ? N : desde + (N / T); // Último hilo toma las filas que sobren

    // Multiplicación clásica de matrices para su rango de filas
    for(int i = desde; i < hasta; i++){
        for(int j = 0; j < N; j++){
            double suma = 0.0;
            for(int k = 0; k < N; k++){
                suma += A[i][k] * B[k][j];
            }
            C[i][j] = suma;
        }
    }

    pthread_exit(NULL); // El hilo finaliza
}

// Inicializa las matrices y les asigna valores
void inicializarMatrices(){
    A = malloc(N * sizeof(double *));
    B = malloc(N * sizeof(double *));
    C = malloc(N * sizeof(double *));

    srand(time(NULL)); // Semilla para generar números aleatorios diferentes en cada ejecución

    for(int i = 0; i < N; i++){
        A[i] = malloc(N * sizeof(double));
        B[i] = malloc(N * sizeof(double));
        C[i] = malloc(N * sizeof(double));

        for(int j = 0; j < N; j++){
            A[i][j] = (double)(rand() % 10);  // Números aleatorios entre 0 y 9
            B[i][j] = (double)(rand() % 10);
            C[i][j] = 0.0; // Resultado inicial en cero
        }
    }
}


// Libera la memoria reservada para las matrices
void liberarMatrices(){
    for(int i = 0; i < N; i++){
        free(A[i]);
        free(B[i]);
        free(C[i]);
    }
    free(A);
    free(B);
    free(C);
}

// Función principal
int main(int argc, char *argv[]){
    // Verifica que se pasen correctamente los parámetros
    if(argc != 3){
        printf("Uso: %s <nro_hilos> <dimension_N>\n", argv[0]);
        return -1;
    }

    T = atoi(argv[1]); // Número de hilos
    N = atoi(argv[2]); // Tamaño de la matriz

    pthread_t threads[T];  // Arreglo para los identificadores de hilos
    int ids[T];            // Arreglo para pasar el ID de cada hilo

    inicializarMatrices(); // Reserva memoria e inicializa A y B

    double inicio = dwalltime(); // Guarda el tiempo de inicio

    // Crea los hilos
    for(int i = 0; i < T; i++){
        ids[i] = i;
        pthread_create(&threads[i], NULL, multiplicar, (void *)&ids[i]);
    }

    // Espera que todos los hilos terminen
    for(int i = 0; i < T; i++){
        pthread_join(threads[i], NULL);
    }

    double fin = dwalltime(); // Tiempo final

    // Muestra el tiempo de ejecución
    printf("Tiempo de ejecución con %d hilos y matriz %d x %d: %.5f segundos\n", T, N, N, fin - inicio);

    liberarMatrices(); // Libera la memoria

    return 0;
}
