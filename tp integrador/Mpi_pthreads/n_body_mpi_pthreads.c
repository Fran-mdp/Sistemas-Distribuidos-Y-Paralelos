// n_body_mpi_pthreads.c
// Versión híbrida MPI + Pthreads para el problema N-cuerpos

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>
#include <mpi.h>

// Constantes físicas y tipos de cuerpos
#define PI 3.141592653589793
#define G 6.673e-11
#define ESTRELLA 0
#define POLVO 1
#define H2 2 // Hidrógeno molecular

// Estructura de cada cuerpo
typedef struct cuerpo {
    float masa;
    float px, py, pz;
    float vx, vy, vz;
    float r, g, b;
    int cuerpo;
} cuerpo_t;

// =====================
// Variables globales (solo para cada proceso)
// =====================
float *fuerza_totalX, *fuerza_totalY, *fuerza_totalZ;
cuerpo_t *cuerpos;          // Subset de cuerpos de este proceso
cuerpo_t *all_cuerpos;      // Todos los cuerpos (para comunicación)
int N, pasos, n_hilos, rank, size;
float delta_tiempo = 1.0f;

// ====================
// Funciones auxiliares
// ====================
double dwalltime() {
    double sec;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

// Inicialización (idéntica a las anteriores)
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
    srand(time(NULL) + start);
    for (cuerpo = 0; cuerpo < local_N; cuerpo++) {
        cuerpos[cuerpo].cuerpo = rand() % 3;
        if (cuerpos[cuerpo].cuerpo == ESTRELLA)
            inicializarEstrella(&cuerpos[cuerpo], start + cuerpo, n);
        else if (cuerpos[cuerpo].cuerpo == POLVO)
            inicializarPolvo(&cuerpos[cuerpo], start + cuerpo, n);
        else
            inicializarH2(&cuerpos[cuerpo], start + cuerpo, n);
    }
    // Opcional: forzar dos cuerpos fijos (solo en rank==0 y rank==1)
    if (rank == 0 && local_N > 0) {
        cuerpos[0].masa = 2.0e2; cuerpos[0].px = 0.0; cuerpos[0].py = 0.0; cuerpos[0].pz = 0.0;
        cuerpos[0].vx = -0.000001; cuerpos[0].vy = -0.000001; cuerpos[0].vz = 0.0;
    }
    if (rank == 1 && local_N > 0) {
        cuerpos[0].masa = 1.0e1; cuerpos[0].px = -1.0; cuerpos[0].py = 0.0; cuerpos[0].pz = 0.0;
        cuerpos[0].vx = 0.0; cuerpos[0].vy = 0.0001; cuerpos[0].vz = 0.0;
    }
}

// ===============================
// Pthreads: paralelismo interno
// ===============================

typedef struct {
    int hilo_id;
    int hilos_totales;
    int N_local;
    int N_global;
    cuerpo_t *cuerpos_local;
    cuerpo_t *cuerpos_all;
    float *fuerza_localX;
    float *fuerza_localY;
    float *fuerza_localZ;
    int start_global; // índice global donde empieza mi bloque
} pthread_args_t;

// Dado el índice plano p, obtiene (cuerpo1, cuerpo2)
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

// Cada hilo calcula las fuerzas entre pares de cuerpos, pero **solo suma fuerza en los cuerpos locales**
void *fuerzas_worker(void *args) {
    pthread_args_t *arg = (pthread_args_t *)args;
    int N_global = arg->N_global;
    int N_local = arg->N_local;
    int hilo_id = arg->hilo_id;
    int hilos_totales = arg->hilos_totales;
    int start_global = arg->start_global;

    int total_pares = (N_global * (N_global - 1)) / 2;
    int pares_por_hilo = total_pares / hilos_totales;
    int inicio = hilo_id * pares_por_hilo;
    int fin = (hilo_id == hilos_totales - 1) ? total_pares : inicio + pares_por_hilo;

    // Inicializar fuerzas locales
    for (int i = 0; i < N_local; i++) {
        arg->fuerza_localX[i] = 0.0;
        arg->fuerza_localY[i] = 0.0;
        arg->fuerza_localZ[i] = 0.0;
    }

    for (int p = inicio; p < fin; p++) {
        int c1, c2;
        get_pair(p, N_global, &c1, &c2);

        // Solo suma fuerzas si el cuerpo pertenece a mi bloque local
        int local_c1 = (c1 >= start_global && c1 < start_global + N_local) ? c1 - start_global : -1;
        int local_c2 = (c2 >= start_global && c2 < start_global + N_local) ? c2 - start_global : -1;

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

        // Solo actualiza los cuerpos locales
        if (local_c1 >= 0) {
            arg->fuerza_localX[local_c1] += dif_X;
            arg->fuerza_localY[local_c1] += dif_Y;
            arg->fuerza_localZ[local_c1] += dif_Z;
        }
        if (local_c2 >= 0) {
            arg->fuerza_localX[local_c2] -= dif_X;
            arg->fuerza_localY[local_c2] -= dif_Y;
            arg->fuerza_localZ[local_c2] -= dif_Z;
        }
    }
    return NULL;
}

// Paraleliza el cálculo de fuerzas usando n_hilos dentro de cada proceso
void calcularFuerzas_pthreads(cuerpo_t *all_cuerpos, cuerpo_t *cuerpos_local, int N_global, int N_local, int n_hilos, float *fuerza_totalX, float *fuerza_totalY, float *fuerza_totalZ, int start_global) {
    pthread_t *hilos = malloc(sizeof(pthread_t) * n_hilos);
    pthread_args_t *args = malloc(sizeof(pthread_args_t) * n_hilos);

    float **fuerzas_locX = malloc(sizeof(float *) * n_hilos);
    float **fuerzas_locY = malloc(sizeof(float *) * n_hilos);
    float **fuerzas_locZ = malloc(sizeof(float *) * n_hilos);

    for (int i = 0; i < n_hilos; i++) {
        fuerzas_locX[i] = calloc(N_local, sizeof(float));
        fuerzas_locY[i] = calloc(N_local, sizeof(float));
        fuerzas_locZ[i] = calloc(N_local, sizeof(float));
    }

    for (int i = 0; i < n_hilos; i++) {
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
    for (int i = 0; i < n_hilos; i++) {
        pthread_join(hilos[i], NULL);
    }

    // Reducción local
    for (int j = 0; j < N_local; j++) {
        fuerza_totalX[j] = 0.0;
        fuerza_totalY[j] = 0.0;
        fuerza_totalZ[j] = 0.0;
        for (int i = 0; i < n_hilos; i++) {
            fuerza_totalX[j] += fuerzas_locX[i][j];
            fuerza_totalY[j] += fuerzas_locY[i][j];
            fuerza_totalZ[j] += fuerzas_locZ[i][j];
        }
    }
    for (int i = 0; i < n_hilos; i++) {
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

// ===================
// Movimiento cuerpos
// ===================
void moverCuerpos(cuerpo_t *cuerpos, int N, int dt, float *fuerza_totalX, float *fuerza_totalY, float *fuerza_totalZ) {
    for (int cuerpo = 0; cuerpo < N; cuerpo++) {
        fuerza_totalX[cuerpo] *= 1 / cuerpos[cuerpo].masa;
        fuerza_totalY[cuerpo] *= 1 / cuerpos[cuerpo].masa;

        cuerpos[cuerpo].vx += fuerza_totalX[cuerpo] * dt;
        cuerpos[cuerpo].vy += fuerza_totalY[cuerpo] * dt;

        cuerpos[cuerpo].px += cuerpos[cuerpo].vx * dt;
        cuerpos[cuerpo].py += cuerpos[cuerpo].vy * dt;

        fuerza_totalX[cuerpo] = 0.0;
        fuerza_totalY[cuerpo] = 0.0;
        fuerza_totalZ[cuerpo] = 0.0;
    }
}

// =================
// MAIN DEL PROGRAMA
// =================
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 5) {
        if (rank == 0)
            printf("Uso: %s <n cuerpos> <DT> <pasos> <n_hilos>\n", argv[0]);
        MPI_Finalize();
        return -1;
    }

    N = atoi(argv[1]);
    delta_tiempo = atof(argv[2]);
    pasos = atoi(argv[3]);
    n_hilos = atoi(argv[4]);

    // División del trabajo entre procesos
    int local_N = N / size;
    int resto = N % size;
    int start_global = rank * local_N + (rank < resto ? rank : resto);
    if (rank < resto)
        local_N += 1; // Los primeros procesos se llevan uno más

    cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * local_N);
    fuerza_totalX = (float *)malloc(sizeof(float) * local_N);
    fuerza_totalY = (float *)malloc(sizeof(float) * local_N);
    fuerza_totalZ = (float *)malloc(sizeof(float) * local_N);

    // Todos los cuerpos globales para la comunicación
    all_cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);

    inicializarCuerpos(cuerpos, N, start_global, local_N);

    // Tipo derivado MPI para la estructura cuerpo_t
    MPI_Datatype MPI_CUERPO;
    MPI_Type_contiguous(sizeof(cuerpo_t)/sizeof(float), MPI_FLOAT, &MPI_CUERPO);
    MPI_Type_commit(&MPI_CUERPO);

    double tIni = 0, tFin = 0;
    if (rank == 0) tIni = dwalltime();

    for (int paso = 0; paso < pasos; paso++) {
        // Todos los procesos comparten sus cuerpos con todos
        MPI_Allgather(cuerpos, local_N * sizeof(cuerpo_t), MPI_BYTE,
                      all_cuerpos, local_N * sizeof(cuerpo_t), MPI_BYTE, MPI_COMM_WORLD);

        // Calcula las fuerzas usando Pthreads
        calcularFuerzas_pthreads(all_cuerpos, cuerpos, N, local_N, n_hilos,
                                fuerza_totalX, fuerza_totalY, fuerza_totalZ, start_global);

        // Mueve los cuerpos locales
        moverCuerpos(cuerpos, local_N, delta_tiempo, fuerza_totalX, fuerza_totalY, fuerza_totalZ);
    }

    if (rank == 0) {
        tFin = dwalltime();
        printf("Tiempo en segundos: %f\n", tFin - tIni);
    }

    free(cuerpos);
    free(fuerza_totalX); free(fuerza_totalY); free(fuerza_totalZ);
    free(all_cuerpos);

    MPI_Type_free(&MPI_CUERPO);
    MPI_Finalize();
    return 0;
}
 