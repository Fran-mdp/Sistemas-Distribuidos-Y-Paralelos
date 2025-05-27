void merge(int izq, int medio, int der) {
    // Calculamos cuántos elementos hay en total entre izq y der
    int n = der - izq + 1;

    // Creamos un arreglo auxiliar para guardar los elementos ordenados temporalmente
    int *aux = malloc(n * sizeof(int));

    // i recorre la primera mitad: desde izq hasta medio
    // j recorre la segunda mitad: desde medio+1 hasta der
    // k es el índice del arreglo aux
    int i = izq;
    int j = medio + 1;
    int k = 0;

    // Mientras haya elementos en ambas mitades...
    while (i <= medio && j <= der) {
        // Comparamos los elementos actuales de la izquierda y la derecha
        // y copiamos el menor al arreglo temporal aux
        if (vector[i] <= vector[j]) {
            aux[k] = vector[i];  // Copiamos el de la izquierda
            i++;                 // Avanzamos en la izquierda
        } else {
            aux[k] = vector[j];  // Copiamos el de la derecha
            j++;                 // Avanzamos en la derecha
        }
        k++; // Avanzamos en el aux
    }

    // Si quedaron elementos sin copiar en la mitad izquierda
    while (i <= medio) {
        aux[k] = vector[i];  // Copiamos el elemento restante
        i++;
        k++;
    }

    // Si quedaron elementos sin copiar en la mitad derecha
    while (j <= der) {
        aux[k] = vector[j];  // Copiamos el elemento restante
        j++;
        k++;
    }

    // Copiamos todos los elementos ordenados del aux al vector original
    for (int p = 0; p < n; p++) {
        vector[izq + p] = aux[p];  // Sobrescribe en el vector original
    }

    // Liberamos la memoria del arreglo auxiliar
    free(aux);
}
