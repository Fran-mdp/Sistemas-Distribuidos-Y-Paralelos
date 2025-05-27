/// n_body_pthreads.c
// Versión paralela usando Pthreads para el problema N-cuerpos (OPTIMIZADA Y PROBADA)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>

// Constantes físicas y tipos de cuerpos
#define PI 3.141592653589793             // Pi para cálculos trigonométricos
#define G 6.673e-11                      // Constante de gravitación universal
#define ESTRELLA 0                       // Identificador para tipo de cuerpo: estrella
#define POLVO 1                          // Identificador para tipo de cuerpo: polvo
#define H2 2                             // Identificador para tipo de cuerpo: hidrógeno molecular

// Estructura que representa un cuerpo
typedef struct cuerpo {
    float masa;          // Masa del cuerpo
    float px, py, pz;    // Posición en X, Y, Z
    float vx, vy, vz;    // Velocidad en X, Y, Z
    float r, g, b;       // Color (para visualización, no afecta cálculo)
    int cuerpo;          // Tipo de cuerpo (ESTRELLA, POLVO, H2)
} cuerpo_t;

// Estructura de argumentos que recibe cada hilo
typedef struct {
    int hilo_id;                // ID del hilo (número de hilo)
    int hilos_totales;          // Total de hilos utilizados
    int N;                      // Número de cuerpos
    cuerpo_t *cuerpos;          // Puntero al array de cuerpos
    float *fuerza_localX;       // Array local de fuerza en X para este hilo
    float *fuerza_localY;       // Array local de fuerza en Y para este hilo
    float *fuerza_localZ;       // Array local de fuerza en Z para este hilo
} pthread_args_t;

// ====================
// Variables globales
// ====================
float *fuerza_totalX, *fuerza_totalY, *fuerza_totalZ;  // Arrays para fuerzas totales sobre cada cuerpo
cuerpo_t *cuerpos;                                     // Array de cuerpos
int N, pasos;                                          // N = número de cuerpos, pasos = pasos de simulación
float delta_tiempo = 1.0f;                             // Paso temporal

// Arrays auxiliares para paralelismo
float **fuerzas_locX;   // Arrays de fuerzas locales en X para cada hilo
float **fuerzas_locY;   // Arrays de fuerzas locales en Y para cada hilo
float **fuerzas_locZ;   // Arrays de fuerzas locales en Z para cada hilo
pthread_args_t *args;   // Array de estructuras de argumentos para hilos
pthread_t *hilos;       // Array de identificadores de hilos
int n_hilos_global;     // Cantidad global de hilos

// ====================
// Funciones auxiliares
// ====================

// Devuelve el tiempo actual en segundos (útil para medir performance)
double dwalltime() {
    double sec;
    struct timeval tv;
    gettimeofday(&tv, NULL);                 // Obtiene el tiempo actual
    sec = tv.tv_sec + tv.tv_usec / 1000000.0; // Convierte a segundos
    return sec;
}

// Inicializa un cuerpo tipo estrella en posición circular
void inicializarEstrella(cuerpo_t *cuerpo, int i, double n) {
    cuerpo->masa = 0.001 * 8;                  // Masa característica de estrella
    cuerpo->px = cos(2 * PI * i / n);          // Posición inicial X sobre círculo
    cuerpo->py = sin(2 * PI * i / n);          // Posición inicial Y sobre círculo
    cuerpo->pz = 0.0;                          // Plano XY (sin Z)
    cuerpo->vx = cuerpo->vy = cuerpo->vz = 0.0;// Velocidades iniciales en 0
    cuerpo->r = 1.0; cuerpo->g = 1.0; cuerpo->b = 1.0; // Color blanco
}
// Inicializa un cuerpo tipo polvo
void inicializarPolvo(cuerpo_t *cuerpo, int i, double n) {
    cuerpo->masa = 0.001 * 4;
    cuerpo->px = cos(2 * PI * i / n);
    cuerpo->py = sin(2 * PI * i / n);
    cuerpo->pz = 0.0;
    cuerpo->vx = cuerpo->vy = cuerpo->vz = 0.0;
    cuerpo->r = 1.0; cuerpo->g = 0.0; cuerpo->b = 0.0; // Color rojo
}
// Inicializa un cuerpo tipo H2 (hidrógeno molecular)
void inicializarH2(cuerpo_t *cuerpo, int i, double n) {
    cuerpo->masa = 0.001;
    cuerpo->px = cos(2 * PI * i / n);
    cuerpo->py = sin(2 * PI * i / n);
    cuerpo->pz = 0.0;
    cuerpo->vx = cuerpo->vy = cuerpo->vz = 0.0;
    cuerpo->r = 1.0; cuerpo->g = 1.0; cuerpo->b = 1.0; // Color blanco
}

// Inicializa todos los cuerpos (N en total)
void inicializarCuerpos(cuerpo_t *cuerpos, int N) {
    int cuerpo;
    double n = N;
    srand(time(NULL)); // Inicializa generador aleatorio
    for (cuerpo = 0; cuerpo < N; cuerpo++) {
        fuerza_totalX[cuerpo] = 0.0;  // Fuerzas totales a 0
        fuerza_totalY[cuerpo] = 0.0;
        fuerza_totalZ[cuerpo] = 0.0;
        cuerpos[cuerpo].cuerpo = rand() % 3; // Tipo aleatorio
        if (cuerpos[cuerpo].cuerpo == ESTRELLA)
            inicializarEstrella(&cuerpos[cuerpo], cuerpo, n);
        else if (cuerpos[cuerpo].cuerpo == POLVO)
            inicializarPolvo(&cuerpos[cuerpo], cuerpo, n);
        else
            inicializarH2(&cuerpos[cuerpo], cuerpo, n);
    }
    // Cuerpos 0 y 1 con parámetros especiales (igual que el código de la cátedra)
    cuerpos[0].masa = 2.0e2; cuerpos[0].px = 0.0; cuerpos[0].py = 0.0; cuerpos[0].pz = 0.0;
    cuerpos[0].vx = -0.000001; cuerpos[0].vy = -0.000001; cuerpos[0].vz = 0.0;
    cuerpos[1].masa = 1.0e1; cuerpos[1].px = -1.0; cuerpos[1].py = 0.0; cuerpos[1].pz = 0.0;
    cuerpos[1].vx = 0.0; cuerpos[1].vy = 0.0001; cuerpos[1].vz = 0.0;
}

// ===============================
// Funciones para paralelización
// ===============================

// Mapea un índice plano p al par único (cuerpo1, cuerpo2) tal que cuerpo1 < cuerpo2
void get_pair(int p, int N, int *c1, int *c2) {
    int sum = 0;
    for (int i = 0; i < N-1; i++) {
        int count = N - i - 1;
        if (p < sum + count) {
            *c1 = i;
            *c2 = i + 1 + (p - sum);
            return;
        }
        sum += count;
    }
}

// Función ejecutada por cada hilo (worker)
void *fuerzas_worker(void *args) {
    pthread_args_t *arg = (pthread_args_t *)args; // Cast a struct de argumentos
    int N = arg->N;
    int hilo_id = arg->hilo_id;
    int hilos_totales = arg->hilos_totales;

    int total_pares = (N * (N - 1)) / 2; // Total de pares de cuerpos
    int pares_por_hilo = total_pares / hilos_totales; // Cuántos pares procesa cada hilo
    int inicio = hilo_id * pares_por_hilo; // Índice inicial de pares para este hilo
    int fin = (hilo_id == hilos_totales - 1) ? total_pares : inicio + pares_por_hilo; // Índice final

    // Poner en cero fuerzas locales de este hilo
    for(int i = 0; i < N; i++) {
        arg->fuerza_localX[i] = 0.0;
        arg->fuerza_localY[i] = 0.0;
        arg->fuerza_localZ[i] = 0.0;
    }

    // Calcular fuerza de cada par asignado a este hilo
    for (int p = inicio; p < fin; p++) {
        int c1, c2;
        get_pair(p, N, &c1, &c2); // Traducir índice plano a par (c1, c2)

        // Si están exactamente en el mismo lugar, se salta (evita división por 0)
        if ((arg->cuerpos[c1].px == arg->cuerpos[c2].px) &&
            (arg->cuerpos[c1].py == arg->cuerpos[c2].py) &&
            (arg->cuerpos[c1].pz == arg->cuerpos[c2].pz))
            continue;

        // Diferencia de posiciones entre cuerpos
        float dif_X = arg->cuerpos[c2].px - arg->cuerpos[c1].px;
        float dif_Y = arg->cuerpos[c2].py - arg->cuerpos[c1].py;
        float dif_Z = arg->cuerpos[c2].pz - arg->cuerpos[c1].pz;

        // Distancia euclídea entre cuerpos
        float distancia = sqrt(dif_X*dif_X + dif_Y*dif_Y + dif_Z*dif_Z);

        // Fuerza gravitatoria entre ambos cuerpos
        float F = (G * arg->cuerpos[c1].masa * arg->cuerpos[c2].masa) / (distancia * distancia);

        // Multiplicar las diferencias por la fuerza (para cada coordenada)
        dif_X *= F;
        dif_Y *= F;
        dif_Z *= F;

        // Acumular fuerza sobre cada cuerpo (acción y reacción)
        arg->fuerza_localX[c1] += dif_X;
        arg->fuerza_localY[c1] += dif_Y;
        arg->fuerza_localZ[c1] += dif_Z;

        arg->fuerza_localX[c2] -= dif_X;
        arg->fuerza_localY[c2] -= dif_Y;
        arg->fuerza_localZ[c2] -= dif_Z;
    }
    return NULL; // Fin del hilo
}

// Función que coordina la paralelización del cálculo de fuerzas
void calcularFuerzas_pthreads(cuerpo_t *cuerpos, int N, int hilos_totales) {
    // Inicializar fuerzas locales a cero para cada hilo
    for(int i=0; i<hilos_totales; i++) {
        for(int j=0; j<N; j++) {
            fuerzas_locX[i][j] = 0.0;
            fuerzas_locY[i][j] = 0.0;
            fuerzas_locZ[i][j] = 0.0;
        }
    }

    // Crear los hilos para que procesen los pares en paralelo
    for(int i=0; i<hilos_totales; i++) {
        args[i].hilo_id = i;
        args[i].hilos_totales = hilos_totales;
        args[i].N = N;
        args[i].cuerpos = cuerpos;
        args[i].fuerza_localX = fuerzas_locX[i];
        args[i].fuerza_localY = fuerzas_locY[i];
        args[i].fuerza_localZ = fuerzas_locZ[i];
        pthread_create(&hilos[i], NULL, fuerzas_worker, &args[i]);
    }

    // Esperar a que todos los hilos terminen su trabajo
    for(int i=0; i<hilos_totales; i++) {
        pthread_join(hilos[i], NULL);
    }

    // Reducción: sumar fuerzas locales de todos los hilos para cada cuerpo
    for(int j=0; j<N; j++) {
        fuerza_totalX[j] = 0.0;
        fuerza_totalY[j] = 0.0;
        fuerza_totalZ[j] = 0.0;
        for(int i=0; i<hilos_totales; i++) {
            fuerza_totalX[j] += fuerzas_locX[i][j];
            fuerza_totalY[j] += fuerzas_locY[i][j];
            fuerza_totalZ[j] += fuerzas_locZ[i][j];
        }
    }
}

// =============
// Mover cuerpos
// =============

// Actualiza la velocidad y la posición de todos los cuerpos, usando las fuerzas calculadas
void moverCuerpos(cuerpo_t *cuerpos, int N, int dt) {
    int cuerpo;
    for (cuerpo = 0; cuerpo < N; cuerpo++) {
        // Calcula aceleración (fuerza / masa) y actualiza velocidad
        fuerza_totalX[cuerpo] *= 1 / cuerpos[cuerpo].masa;
        fuerza_totalY[cuerpo] *= 1 / cuerpos[cuerpo].masa;

        cuerpos[cuerpo].vx += fuerza_totalX[cuerpo] * dt; // Nueva velocidad X
        cuerpos[cuerpo].vy += fuerza_totalY[cuerpo] * dt; // Nueva velocidad Y

        cuerpos[cuerpo].px += cuerpos[cuerpo].vx * dt;    // Nueva posición X
        cuerpos[cuerpo].py += cuerpos[cuerpo].vy * dt;    // Nueva posición Y

        // Limpia fuerzas para la próxima iteración
        fuerza_totalX[cuerpo] = 0.0;
        fuerza_totalY[cuerpo] = 0.0;
        fuerza_totalZ[cuerpo] = 0.0;
    }
}

// =========================
// Liberación de memoria
// =========================
void liberarMemoria() {
    free(cuerpos);
    free(fuerza_totalX); free(fuerza_totalY); free(fuerza_totalZ);
    for(int i=0; i<n_hilos_global; i++) {
        free(fuerzas_locX[i]);
        free(fuerzas_locY[i]);
        free(fuerzas_locZ[i]);
    }
    free(fuerzas_locX); free(fuerzas_locY); free(fuerzas_locZ);
    free(args); free(hilos);
}

// =================
// MAIN DEL PROGRAMA
// =================

int main(int argc, char *argv[]) {
    if (argc < 5) {
        printf("Uso: %s <n cuerpos> <DT> <pasos> <n_hilos>\n", argv[0]);
        return -1;
    }

    N = atoi(argv[1]);                     // Número de cuerpos
    delta_tiempo = atof(argv[2]);          // Paso temporal
    pasos = atoi(argv[3]);                 // Cantidad de pasos de simulación
    int n_hilos = atoi(argv[4]);           // Cantidad de hilos (cores)
    n_hilos_global = n_hilos;

    // Reserva memoria para cuerpos y fuerzas
    cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);
    fuerza_totalX = (float *)malloc(sizeof(float) * N);
    fuerza_totalY = (float *)malloc(sizeof(float) * N);
    fuerza_totalZ = (float *)malloc(sizeof(float) * N);

    // Reserva memoria para las fuerzas locales y estructuras de hilos
    fuerzas_locX = malloc(sizeof(float*) * n_hilos_global);
    fuerzas_locY = malloc(sizeof(float*) * n_hilos_global);
    fuerzas_locZ = malloc(sizeof(float*) * n_hilos_global);
    args        = malloc(sizeof(pthread_args_t) * n_hilos_global);
    hilos       = malloc(sizeof(pthread_t) * n_hilos_global);

    // Reserva espacio para los arrays de fuerzas de cada hilo
    for(int i=0; i<n_hilos_global; i++) {
        fuerzas_locX[i] = calloc(N, sizeof(float));
        fuerzas_locY[i] = calloc(N, sizeof(float));
        fuerzas_locZ[i] = calloc(N, sizeof(float));
    }

    inicializarCuerpos(cuerpos, N); // Inicializa posiciones y masas de los cuerpos

    double tIni = dwalltime(); // Marca el tiempo de inicio

    // Bucle principal: para cada paso de simulación
    for (int paso = 0; paso < pasos; paso++) {
        calcularFuerzas_pthreads(cuerpos, N, n_hilos);     // Calcula todas las fuerzas (en paralelo)
        moverCuerpos(cuerpos, N, delta_tiempo);            // Actualiza posiciones y velocidades
    }

    double tFin = dwalltime(); // Marca el tiempo de fin
    printf("Tiempo en segundos: %f\n", tFin - tIni); // Muestra el tiempo de simulación

    liberarMemoria(); // Libera toda la memoria antes de salir

    return 0; // Fin
}
