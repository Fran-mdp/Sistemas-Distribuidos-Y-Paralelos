// n_body_simple.c
// Versión secuencial del problema de N-cuerpos

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

// Constantes
#define PI 3.141592653589793
#define G 6.673e-11
#define ESTRELLA 0
#define POLVO 1
#define H2 2 // Hidrógeno molecular

// Estructura de cada cuerpo
typedef struct cuerpo {
    float masa;
    float px, py, pz;     // posición
    float vx, vy, vz;     // velocidad
    float r, g, b;        // color (no se usa en la versión NO-GL)
    int cuerpo;           // tipo de cuerpo
} cuerpo_t;

// Variables globales
float *fuerza_totalX, *fuerza_totalY, *fuerza_totalZ;
cuerpo_t *cuerpos;

// Parámetros de simulación
int delta_tiempo = 1.0f; // Intervalo de tiempo por paso
int pasos;
int N;

// ====================
// Funciones auxiliares
// ====================

// Devuelve el tiempo actual en segundos (para medir performance)
double dwalltime() {
    double sec;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

// Inicializa un cuerpo como "estrella"
void inicializarEstrella(cuerpo_t *cuerpo, int i, double n) {
    cuerpo->masa = 0.001*8;
    // Ubicación simple en un círculo para ejemplo
    cuerpo->px = cos(2*PI*i/n);
    cuerpo->py = sin(2*PI*i/n);
    cuerpo->pz = 0.0;
    cuerpo->vx = cuerpo->vy = cuerpo->vz = 0.0;
    cuerpo->r = 1.0; cuerpo->g = 1.0; cuerpo->b = 1.0;
}

// Inicializa un cuerpo como "polvo"
void inicializarPolvo(cuerpo_t *cuerpo, int i, double n) {
    cuerpo->masa = 0.001*4;
    cuerpo->px = cos(2*PI*i/n);
    cuerpo->py = sin(2*PI*i/n);
    cuerpo->pz = 0.0;
    cuerpo->vx = cuerpo->vy = cuerpo->vz = 0.0;
    cuerpo->r = 1.0; cuerpo->g = 0.0; cuerpo->b = 0.0;
}

// Inicializa un cuerpo como "H2"
void inicializarH2(cuerpo_t *cuerpo, int i, double n) {
    cuerpo->masa = 0.001;
    cuerpo->px = cos(2*PI*i/n);
    cuerpo->py = sin(2*PI*i/n);
    cuerpo->pz = 0.0;
    cuerpo->vx = cuerpo->vy = cuerpo->vz = 0.0;
    cuerpo->r = 1.0; cuerpo->g = 1.0; cuerpo->b = 1.0;
}

// Inicializa todos los cuerpos y fuerzas
void inicializarCuerpos(cuerpo_t *cuerpos, int N) {
    int cuerpo;
    double n = N;
    srand(time(NULL));
    for (cuerpo = 0; cuerpo < N; cuerpo++) {
        fuerza_totalX[cuerpo] = 0.0;
        fuerza_totalY[cuerpo] = 0.0;
        fuerza_totalZ[cuerpo] = 0.0;
        cuerpos[cuerpo].cuerpo = rand() % 3;
        if (cuerpos[cuerpo].cuerpo == ESTRELLA)
            inicializarEstrella(&cuerpos[cuerpo], cuerpo, n);
        else if (cuerpos[cuerpo].cuerpo == POLVO)
            inicializarPolvo(&cuerpos[cuerpo], cuerpo, n);
        else
            inicializarH2(&cuerpos[cuerpo], cuerpo, n);
    }

    // Asigna manualmente los dos primeros cuerpos como ejemplo
    cuerpos[0].masa = 2.0e2; cuerpos[0].px = 0.0; cuerpos[0].py = 0.0; cuerpos[0].pz = 0.0;
    cuerpos[0].vx = -0.000001; cuerpos[0].vy = -0.000001; cuerpos[0].vz = 0.0;
    cuerpos[1].masa = 1.0e1; cuerpos[1].px = -1.0; cuerpos[1].py = 0.0; cuerpos[1].pz = 0.0;
    cuerpos[1].vx = 0.0; cuerpos[1].vy = 0.0001; cuerpos[1].vz = 0.0;
}

// =====================
// Simulación secuencial
// =====================

// Calcula todas las fuerzas gravitatorias (doble for secuencial)
void calcularFuerzas(cuerpo_t *cuerpos, int N, int dt) {
    int cuerpo1, cuerpo2;
    float dif_X, dif_Y, dif_Z, distancia, F;

    for (cuerpo1 = 0; cuerpo1 < N - 1; cuerpo1++) {
        for (cuerpo2 = cuerpo1 + 1; cuerpo2 < N; cuerpo2++) {
            if ((cuerpos[cuerpo1].px == cuerpos[cuerpo2].px) &&
                (cuerpos[cuerpo1].py == cuerpos[cuerpo2].py) &&
                (cuerpos[cuerpo1].pz == cuerpos[cuerpo2].pz))
                continue;

            dif_X = cuerpos[cuerpo2].px - cuerpos[cuerpo1].px;
            dif_Y = cuerpos[cuerpo2].py - cuerpos[cuerpo1].py;
            dif_Z = cuerpos[cuerpo2].pz - cuerpos[cuerpo1].pz;

            distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);

            F = (G * cuerpos[cuerpo1].masa * cuerpos[cuerpo2].masa) / (distancia * distancia);

            dif_X *= F;
            dif_Y *= F;
            dif_Z *= F;

            fuerza_totalX[cuerpo1] += dif_X;
            fuerza_totalY[cuerpo1] += dif_Y;
            fuerza_totalZ[cuerpo1] += dif_Z;

            fuerza_totalX[cuerpo2] -= dif_X;
            fuerza_totalY[cuerpo2] -= dif_Y;
            fuerza_totalZ[cuerpo2] -= dif_Z;
        }
    }
}

// Mueve todos los cuerpos en función de las fuerzas
void moverCuerpos(cuerpo_t *cuerpos, int N, int dt) {
    int cuerpo;
    for (cuerpo = 0; cuerpo < N; cuerpo++) {
        fuerza_totalX[cuerpo] *= 1 / cuerpos[cuerpo].masa;
        fuerza_totalY[cuerpo] *= 1 / cuerpos[cuerpo].masa;
        // fuerza_totalZ[cuerpo] *= 1 / cuerpos[cuerpo].masa;

        cuerpos[cuerpo].vx += fuerza_totalX[cuerpo] * dt;
        cuerpos[cuerpo].vy += fuerza_totalY[cuerpo] * dt;
        // cuerpos[cuerpo].vz += fuerza_totalZ[cuerpo] * dt;

        cuerpos[cuerpo].px += cuerpos[cuerpo].vx * dt;
        cuerpos[cuerpo].py += cuerpos[cuerpo].vy * dt;
        // cuerpos[cuerpo].pz += cuerpos[cuerpo].vz * dt;

        fuerza_totalX[cuerpo] = 0.0;
        fuerza_totalY[cuerpo] = 0.0;
        fuerza_totalZ[cuerpo] = 0.0;
    }
}

// =================
// MAIN DEL PROGRAMA
// =================

int main(int argc, char *argv[]) {
    if (argc < 4) {
        printf("Uso: %s <n cuerpos> <DT> <pasos>\n", argv[0]);
        return -1;
    }

    N = atoi(argv[1]);
    delta_tiempo = atof(argv[2]);
    pasos = atoi(argv[3]);

    // Reservar memoria
    cuerpos = (cuerpo_t *)malloc(sizeof(cuerpo_t) * N);
    fuerza_totalX = (float *)malloc(sizeof(float) * N);
    fuerza_totalY = (float *)malloc(sizeof(float) * N);
    fuerza_totalZ = (float *)malloc(sizeof(float) * N);

    // Inicializar
    inicializarCuerpos(cuerpos, N);

    double tIni = dwalltime();

    // Simulación principal
    for (int paso = 0; paso < pasos; paso++) {
        calcularFuerzas(cuerpos, N, delta_tiempo);
        moverCuerpos(cuerpos, N, delta_tiempo);
    }

    double tFin = dwalltime();
    printf("Tiempo en segundos: %f\n", tFin - tIni);

    // Liberar memoria
    free(cuerpos);
    free(fuerza_totalX);
    free(fuerza_totalY);
    free(fuerza_totalZ);

    return 0;
}
