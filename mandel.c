#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>

#define DEBUG 1

#define X_RESN 256 // Resolución en el eje X
#define Y_RESN 256 // Resolución en el eje Y

#define X_MIN -2.0
#define X_MAX 2.0
#define Y_MIN -2.0
#define Y_MAX 2.0

#define maxIterations 1000 //máximo número de iteraciones para cada punto en el conjunto de Mandelbrot

//definición de una estructura para representar números complejos (real e imaginario)
typedef struct complextype {
    float real, imag; //parte real e imaginaria del número complejo
} Compl;

//función auxiliar para calcular el tiempo en segundos entre dos eventos usando struct timeval
static inline double get_seconds(struct timeval t_ini, struct timeval t_end) {
    return (t_end.tv_usec - t_ini.tv_usec) / 1E6 + (t_end.tv_sec - t_ini.tv_sec); // Calcular diferencia de tiempo en segundos
}

int main(int argc, char **argv) {
    int rank, num_procs;  //rank: Identificador único del proceso, num_procs: número total de procesos
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //obtener el identificador (rank) del proceso actual
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); //obtener el número total de procesos

    struct timeval t_total1, t_total2, t_comm1, t_comm2, t_comp1, t_comp2; //variables para medir los tiempos

    int *vres = NULL, **res = NULL;  //variables para almacenar los resultados finales
    int *resultado_local = NULL; //resultado local para cada proceso
    int *array_filas_por_proceso = NULL, *inicio_filas = NULL; //arrays para manejar las filas distribuidas
    int *indice_filas_locales = NULL; //índices de filas locales
    int *filas_totales = NULL; //filas totales

    int resto = Y_RESN % num_procs; //resto de filas cuando Y_RESN no es divisible por num_procs
    int filas_por_proceso = Y_RESN / num_procs; //número de filas que le toca a cada proceso

    //asignación dinámica de memoria para las estructuras de filas
    array_filas_por_proceso = (int *)malloc(num_procs * sizeof(int));
    inicio_filas = (int *)malloc(num_procs * sizeof(int));

    int desplazamiento = 0; //variable para llevar el conteo de filas procesadas
    for (int i = 0; i < num_procs; i++) {
        array_filas_por_proceso[i] = filas_por_proceso + (i < resto ? 1 : 0); //se determina el número de filas por proceso
        inicio_filas[i] = desplazamiento; //se guarda la fila de inicio para cada proceso
        desplazamiento += array_filas_por_proceso[i]; //se actualiza el desplazamiento
    }

    if (rank == 0) { //
        filas_totales = (int *)malloc(Y_RESN * sizeof(int));
        for (int i = 0; i < Y_RESN; i++) {
            filas_totales[i] = i; //se asignan los índices de las filas
        }
    }

    // Calcular cuántas filas le corresponden al proceso actual
    int num_filas_locales = array_filas_por_proceso[rank];
    indice_filas_locales = (int *)malloc(num_filas_locales * sizeof(int));
    resultado_local = (int *)malloc(num_filas_locales * X_RESN * sizeof(int));

    gettimeofday(&t_total1, NULL); //se inicia el contador de tiempo total


    gettimeofday(&t_comm1, NULL); //se iniciar el contador de tiempo de comunicación
    MPI_Scatterv(filas_totales, array_filas_por_proceso, inicio_filas, MPI_INT,
                 indice_filas_locales, num_filas_locales, MPI_INT,
                 0, MPI_COMM_WORLD); //se distribuyen las filas entre los procesos
    gettimeofday(&t_comm2, NULL); //fin del contador de comunicación


    gettimeofday(&t_comp1, NULL); //se inicia el contador de tiempo de computación
    for (int r = 0; r < num_filas_locales; r++) {
        int i = indice_filas_locales[r]; //se obtiene el índice de fila actual
        for (int j = 0; j < X_RESN; j++) { // Para cada columna en la fila
            Compl z = {0.0, 0.0};
            Compl c = {
                X_MIN + j * (X_MAX - X_MIN) / X_RESN,
                Y_MAX - i * (Y_MAX - Y_MIN) / Y_RESN
            };
            int k = 0; //contador de iteraciones
            float lengthsq, temp;

            //se itera hasta alcanzar el límite de escape o el máximo de iteraciones
            do {
                temp = z.real * z.real - z.imag * z.imag + c.real;
                z.imag = 2.0 * z.real * z.imag + c.imag;
                z.real = temp;
                lengthsq = z.real * z.real + z.imag * z.imag;
                k++;
            } while (lengthsq < 4.0 && k < maxIterations);

            resultado_local[r * X_RESN + j] = (k >= maxIterations) ? 0 : k; //se almacena el resultado de la iteración
        }
    }
    gettimeofday(&t_comp2, NULL); //fin del contador de tiempo de computación

    if (rank == 0) { //solo el proceso 0 inicializa la memoria para almacenar los resultados finales
        vres = (int *)malloc(Y_RESN * X_RESN * sizeof(int));
        res = (int **)malloc(Y_RESN * sizeof(int *));
        for (int i = 0; i < Y_RESN; i++) {
            res[i] = vres + i * X_RESN;
        }
    }

    int *filas_recibidas = (int *)malloc(num_procs * sizeof(int));
    int *inicio_resultado = (int *)malloc(num_procs * sizeof(int));
    for (int i = 0; i < num_procs; i++) {
        filas_recibidas[i] = array_filas_por_proceso[i] * X_RESN; //se calcula el tamaño de los datos a recibir por proceso
        inicio_resultado[i] = inicio_filas[i] * X_RESN; //se calcula el desplazamiento en los resultados
    }

    gettimeofday(&t_comm1, NULL); //se inicia el contador de comunicación
    MPI_Gatherv(resultado_local, num_filas_locales * X_RESN, MPI_INT,
                vres, filas_recibidas, inicio_resultado, MPI_INT,
                0, MPI_COMM_WORLD);
    gettimeofday(&t_comm2, NULL); //fin del contador de comunicación

    gettimeofday(&t_total2, NULL); //fin del contador de tiempo total

    double t_total = get_seconds(t_total1, t_total2); //tiempo total
    double t_comm = get_seconds(t_comm1, t_comm2); //tiempo de comunicación
    double t_comp = get_seconds(t_comp1, t_comp2); //tiempo de computación

    //cada proceso imprime su tiempo de computación y comunicación
    fprintf(stderr, "RANK %d:\n  Tiempo computación = %.6lf s\n  Tiempo comunicación = %.6lf s\n",
            rank, t_comp, t_comm);

    //solo el proceso 0 imprime el tiempo total
    if (rank == 0) {
        fprintf(stderr, "(PERF) Time (seconds) = %lf\n", t_total);
    }

    //el proceso 0 imprime el resultado final
    if (rank == 0 && DEBUG) {
        for (int i = 0; i < Y_RESN; i++) {
            for (int j = 0; j < X_RESN; j++) {
                printf("%3d ", res[i][j]); //se imprime cada valor del conjunto de Mandelbrot
            }
            printf("\n");
        }
    }

    //liberacióm de memoria
    if (rank == 0) {
        free(vres);
        free(res);
        free(filas_totales);
    }


    free(resultado_local);
    free(indice_filas_locales);
    free(array_filas_por_proceso);
    free(inicio_filas);
    free(filas_recibidas);
    free(inicio_resultado);

    MPI_Finalize();
    return 0;
}
