/// n_body_mpi_pthreads.c
// Versión híbrida MPI + Pthreads para el problema N-cuerpos

#include <stdio.h>           // Entrada/salida estándar
#include <stdlib.h>          // Para malloc, free, atoi, etc
#include <math.h>            // Funciones matemáticas como sqrt, cos, etc
#include <time.h>            // Para la semilla de rand()
#include <sys/time.h>        // Para medir el tiempo de ejecución
#include <pthread.h>         // Pthreads (hilos)
#include <mpi.h>             // MPI

#define PI 3.141592653589793
#define G 6.673e-11           // Constante gravitacional
#define ESTRELLA 0
#define POLVO 1
#define H2 2                  // Tipos de cuerpos

// Estructura que representa un cuerpo en el sistema
typedef struct cuerpo {
    float masa;              // Masa
    float px, py, pz;        // Posición en x, y, z
    float vx, vy, vz;        // Velocidad en x, y, z
    float r, g, b;           // Color (no se usa en cálculo)
    int cuerpo;              // Tipo de cuerpo (ESTRELLA, POLVO o H2)
} cuerpo_t;

// Variables globales principales para cada proceso
float *fuerza_totalX, *fuerza_totalY, *fuerza_totalZ;  // Fuerzas sobre cada cuerpo local
cuerpo_t *cuerpos;                                     // Vector de cuerpos locales
cuerpo_t *all_cuerpos;                                 // Vector con todos los cuerpos (de todos los procesos)
int N, pasos, n_hilos, rank, size;                     // Parámetros generales y de MPI
float delta_tiempo = 1.0f;                             // Paso temporal (puede cambiarse)

// --------------- Función auxiliar para medir el tiempo ---------------
double dwalltime() {
    double sec;
    struct timeval tv;
    gettimeofday(&tv, NULL);                // Obtiene el tiempo actual
    sec = tv.tv_sec + tv.tv_usec / 1000000.0; // Lo pasa a segundos con decimales
    return sec;
}

// -------------------- Inicialización de cuerpos -----------------------
void inicializarEstrella(cuerpo_t *cuerpo, int i, double n) {
    cuerpo->masa = 0.001 * 8;
    cuerpo->px = cos(2 * PI * i / n);
    cuerpo->py = sin(2 * PI * i / n);
    cuerpo->pz = 0.0;
    cuerpo->vx = cuerpo->vy = cuerpo->vz = 0.0;
    cuerpo->r = 1.0; cuerpo->g = 1.0; cuerpo->b = 1.0;
}
void inicializarPolvo(cuerpo_t *cuerpo, int i, double n) {
    cuerpo->masa = 0.001 * 4;
    cuerpo->px = cos(2 * PI * i / n);
    cuerpo->py = sin(2 * PI * i / n);
    cuerpo->pz = 0.0;
    cuerpo->vx = cuerpo->vy = cuerpo->vz = 0.0;
    cuerpo->r = 1.0; cuerpo->g = 0.0; cuerpo->b = 0.0;
}
void inicializarH2(cuerpo_t *cuerpo, int i, double n) {
    cuerpo->masa = 0.001;
    cuerpo->px = cos(2 * PI * i / n);
    cuerpo->py = sin(2 * PI * i / n);
    cuerpo->pz = 0.0;
    cuerpo->vx = cuerpo->vy = cuerpo->vz = 0.0;
    cuerpo->r = 1.0; cuerpo->g = 1.0; cuerpo->b = 1.0;
}
void inicializarCuerpos(cuerpo_t *cuerpos, int N_global, int start, int local_N) {
    int cuerpo;
    double n = N_global;
    srand(time(NULL) + start);                  // Semilla distinta por proceso
    for (cuerpo = 0; cuerpo < local_N; cuerpo++) {
        cuerpos[cuerpo].cuerpo = rand() % 3;    // Tipo aleatorio
        if (cuerpos[cuerpo].cuerpo == ESTRELLA)
            inicializarEstrella(&cuerpos[cuerpo], start + cuerpo, n);
        else if (cuerpos[cuerpo].cuerpo == POLVO)
            inicializarPolvo(&cuerpos[cuerpo], start + cuerpo, n);
        else
            inicializarH2(&cuerpos[cuerpo], start + cuerpo, n);
    }
    // Inicialización de dos cuerpos especiales (solo para rank 0 y 1)
    if (rank == 0 && local_N > 0) {
        cuerpos[0].masa = 2.0e2; cuerpos[0].px = 0.0; cuerpos[0].py = 0.0; cuerpos[0].pz = 0.0;
        cuerpos[0].vx = -0.000001; cuerpos[0].vy = -0.000001; cuerpos[0].vz = 0.0;
    }
    if (rank == 1 && local_N > 0) {
        cuerpos[0].masa = 1.0e1; cuerpos[0].px = -1.0; cuerpos[0].py = 0.0; cuerpos[0].pz = 0.0;
        cuerpos[0].vx = 0.0; cuerpos[0].vy = 0.0001; cuerpos[0].vz = 0.0;
    }
}

// ------------------------ Pthreads: argumentos ------------------------
typedef struct {
    int hilo_id;                // ID de hilo (0, 1, ..., n_hilos-1)
    int hilos_totales;          // Cantidad de hilos
    int N_local;                // Cuerpos locales
    int N_global;               // Cuerpos totales
    cuerpo_t *cuerpos_local;    // Vector de cuerpos locales
    cuerpo_t *cuerpos_all;      // Vector de todos los cuerpos (global)
    float *fuerza_localX;
    float *fuerza_localY;
    float *fuerza_localZ;
    int start_global;           // Índice global donde arranca mi bloque local
} pthread_args_t;

// ------------- Mapeo de índice plano de pares a (c1,c2) únicos --------
void get_pair(int p, int N, int *c1, int *c2) {
    int sum = 0;
    for (int i = 0; i < N - 1; i++) {
        int count = N - i - 1;
        if (p < sum + count) {
            *c1 = i;
            *c2 = i + 1 + (p - sum);
            return;
        }
        sum += count;
    }
}

// --------- Función que ejecuta cada hilo: calcula fuerzas de pares asignados ---------
void *fuerzas_worker(void *args) {
    pthread_args_t *arg = (pthread_args_t *)args;
    int N_global = arg->N_global;
    int N_local = arg->N_local;
    int hilo_id = arg->hilo_id;
    int hilos_totales = arg->hilos_totales;
    int start_global = arg->start_global;

    int total_pares = (N_global * (N_global - 1)) / 2;         // Total de pares posibles
    int pares_por_hilo = total_pares / hilos_totales;          // Cantidad de pares asignados por hilo
    int inicio = hilo_id * pares_por_hilo;                     // Índice inicial de pares
    int fin = (hilo_id == hilos_totales - 1) ? total_pares : inicio + pares_por_hilo; // Índice final

    for (int i = 0; i < N_local; i++) {                        // Inicializa fuerzas locales a cero
        arg->fuerza_localX[i] = 0.0;
        arg->fuerza_localY[i] = 0.0;
        arg->fuerza_localZ[i] = 0.0;
    }

    for (int p = inicio; p < fin; p++) {                       // Para cada par asignado
        int c1, c2;
        get_pair(p, N_global, &c1, &c2);                       // Obtiene qué cuerpos corresponden al par
        int local_c1 = (c1 >= start_global && c1 < start_global + N_local) ? c1 - start_global : -1; // Si c1 está en este proceso
        int local_c2 = (c2 >= start_global && c2 < start_global + N_local) ? c2 - start_global : -1; // Si c2 está en este proceso

        // Si los dos cuerpos están en la misma posición, ignora el par
        if ((arg->cuerpos_all[c1].px == arg->cuerpos_all[c2].px) &&
            (arg->cuerpos_all[c1].py == arg->cuerpos_all[c2].py) &&
            (arg->cuerpos_all[c1].pz == arg->cuerpos_all[c2].pz))
            continue;

        float dif_X = arg->cuerpos_all[c2].px - arg->cuerpos_all[c1].px;
        float dif_Y = arg->cuerpos_all[c2].py - arg->cuerpos_all[c1].py;
        float dif_Z = arg->cuerpos_all[c2].pz - arg->cuerpos_all[c1].pz;
        float distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);
        float F = (G * arg->cuerpos_all[c1].masa * arg->cuerpos_all[c2].masa) / (distancia * distancia);

        dif_X *= F;
        dif_Y *= F;
        dif_Z *= F;

        if (local_c1 >= 0) {                      // Si c1 pertenece a este proceso, suma fuerza a ese cuerpo
            arg->fuerza_localX[local_c1] += dif_X;
            arg->fuerza_localY[local_c1] += dif_Y;
            arg->fuerza_localZ[local_c1] += dif_Z;
        }
        if (local_c2 >= 0) {                      // Si c2 pertenece a este proceso, resta fuerza a ese cuerpo
            arg->fuerza_localX[local_c2] -= dif_X;
            arg->fuerza_localY[local_c2] -= dif_Y;
            arg->fuerza_localZ[local_c2] -= dif_Z;
        }
    }
    return NULL;
}

// --------- Lanza los hilos y suma/reduce resultados en fuerzas globales ---------
void calcularFuerzas_pthreads(
    cuerpo_t *all_cuerpos, cuerpo_t *cuerpos_local,
    int N_global, int N_local, int n_hilos,
    float *fuerza_totalX, float *fuerza_totalY, float *fuerza_totalZ,
    int start_global
) {
    pthread_t *hilos = malloc(sizeof(pthread_t) * n_hilos);
    pthread_args_t *args = malloc(sizeof(pthread_args_t) * n_hilos);
    float **fuerzas_locX = malloc(sizeof(float *) * n_hilos);
    float **fuerzas_locY = malloc(sizeof(float *) * n_hilos);
    float **fuerzas_locZ = malloc(sizeof(float *) * n_hilos);

    for (int i = 0; i < n_hilos; i++) {                              // Crea arrays para cada hilo
        fuerzas_locX[i] = calloc(N_local, sizeof(float));
        fuerzas_locY[i] = calloc(N_local, sizeof(float));
        fuerzas_locZ[i] = calloc(N_local, sizeof(float));
    }

    for (int i = 0; i < n_hilos; i++) {                              // Inicializa argumentos y lanza hilos
        args[i].hilo_id = i;
        args[i].hilos_totales = n_hilos;
        args[i].N_global = N_global;
        args[i].N_local = N_local;
        args[i].start_global = start_global;
        args[i].cuerpos_local = cuerpos_local;
        args[i].cuerpos_all = all_cuerpos;
        args[i].fuerza_localX = fuerzas_locX[i];
        args[i].fuerza_localY = fuerzas_locY[i];
        args[i].fuerza_localZ = fuerzas_locZ[i];
        pthread_create(&hilos[i], NULL, fuerzas_worker, &args[i]);
    }
    for (int i = 0; i < n_hilos; i++) {                              // Espera a que terminen
        pthread_join(hilos[i], NULL);
    }
    for (int j = 0; j < N_local; j++) {                              // Reducción: suma resultados de todos los hilos
        fuerza_totalX[j] = 0.0;
        fuerza_totalY[j] = 0.0;
        fuerza_totalZ[j] = 0.0;
        for (int i = 0; i < n_hilos; i++) {
            fuerza_totalX[j] += fuerzas_locX[i][j];
            fuerza_totalY[j] += fuerzas_locY[i][j];
            fuerza_totalZ[j] += fuerzas_locZ[i][j];
        }
    }
    for (int i = 0; i < n_hilos; i++) {                              // Libera memoria auxiliar de cada hilo
        free(fuerzas_locX[i]);
        free(fuerzas_locY[i]);
        free(fuerzas_locZ[i]);
    }
    free(fuerzas_locX);
    free(fuerzas_locY);
    free(fuerzas_locZ);
    free(hilos);
    free(args);
}

// --------- Movimiento de los cuerpos locales ---------
void moverCuerpos(
    cuerpo_t *cuerpos, int N, int dt,
    float *fuerza_totalX, float *fuerza_totalY, float *fuerza_totalZ
) {
    for (int cuerpo = 0; cuerpo < N; cuerpo++) {
        fuerza_totalX[cuerpo] *= 1 / cuerpos[cuerpo].masa; // Acelera según masa
        fuerza_totalY[cuerpo] *= 1 / cuerpos[cuerpo].masa;

        cuerpos[cuerpo].vx += fuerza_totalX[cuerpo] * dt;  // Actualiza velocidad
        cuerpos[cuerpo].vy += fuerza_totalY[cuerpo] * dt;

        cuerpos[cuerpo].px += cuerpos[cuerpo].vx * dt;     // Actualiza posición
        cuerpos[cuerpo].py += cuerpos[cuerpo].vy * dt;

        fuerza_totalX[cuerpo] = 0.0;                       // Limpia fuerzas para el siguiente paso
        fuerza_totalY[cuerpo] = 0.0;
        fuerza_totalZ[cuerpo] = 0.0;
    }
}

// --------- Libera la memoria reservada ---------
void liberar_memoria() {
    free(cuerpos);
    free(fuerza_totalX);
    free(fuerza_totalY);
    free(fuerza_totalZ);
    free(all_cuerpos);
}

// ================== MAIN DEL PROGRAMA ==================
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);                        // Inicializa el entorno MPI

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);          // rank = número de proceso (0 a N-1)
    MPI_Comm_size(MPI_COMM_WORLD, &size);          // size = cantidad total de procesos

    if (argc < 5) {                                // Chequea argumentos
        if (rank == 0)
            printf("Uso: %s <n cuerpos> <DT> <pasos> <n_hilos>\n", argv[0]);
        MPI_Finalize();
        return -1;
    }

    N = atoi(argv[1]);                             // Cantidad total de cuerpos
    delta_tiempo = atof(argv[2]);                  // Paso de tiempo de simulación
    pasos = atoi(argv[3]);                         // Cantidad de pasos
    n_hilos = atoi(argv[4]);                       // Hilos pthreads por proceso

    // Divide los cuerpos de manera balanceada entre los procesos
    int local_N = N / size;                        // Cuerpos mínimos por proceso
    int resto = N % size;                          // Los primeros 'resto' procesos tienen 1 cuerpo extra
    int start_global = rank * local_N + (rank < resto ? rank : resto); // Índice global de arranque
    if (rank < resto)
        local_N += 1;

    // Reserva memoria para los arrays principales
    cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * local_N);
    fuerza_totalX = (float *)malloc(sizeof(float) * local_N);
    fuerza_totalY = (float *)malloc(sizeof(float) * local_N);
    fuerza_totalZ = (float *)malloc(sizeof(float) * local_N);
    all_cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);

    inicializarCuerpos(cuerpos, N, start_global, local_N); // Inicializa los cuerpos de este proceso

    double tIni = 0, tFin = 0;                   // Variables para medir tiempo de ejecución
    if (rank == 0) tIni = dwalltime();

    for (int paso = 0; paso < pasos; paso++) {   // Bucle principal de la simulación
        // Sincroniza todos los cuerpos de todos los procesos (copia todos en all_cuerpos)
        MPI_Allgather(cuerpos, local_N * sizeof(cuerpo_t), MPI_BYTE,
                      all_cuerpos, local_N * sizeof(cuerpo_t), MPI_BYTE, MPI_COMM_WORLD);

        // Calcula las fuerzas entre todos los cuerpos (en paralelo con hilos)
        calcularFuerzas_pthreads(all_cuerpos, cuerpos, N, local_N, n_hilos,
                                fuerza_totalX, fuerza_totalY, fuerza_totalZ, start_global);

        // Mueve los cuerpos locales (integración de ecuaciones de movimiento)
        moverCuerpos(cuerpos, local_N, delta_tiempo, fuerza_totalX, fuerza_totalY, fuerza_totalZ);
    }

    if (rank == 0) {                        // Solo el proceso maestro imprime el tiempo total
        tFin = dwalltime();
        printf("Tiempo en segundos: %f\n", tFin - tIni);
    }

    liberar_memoria();                      // Libera memoria reservada

    MPI_Finalize();                         // Finaliza el entorno MPI
    return 0;
}
