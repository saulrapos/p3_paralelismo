#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>

#define DEBUG 1

#define X_RESN 1024
#define Y_RESN 1024

#define X_MIN -2.0
#define X_MAX 2.0
#define Y_MIN -2.0
#define Y_MAX 2.0

#define maxIterations 1000

typedef struct complextype {
    float real, imag;
} Compl;

static inline double get_seconds(struct timeval t_ini, struct timeval t_end) {
    return (t_end.tv_usec - t_ini.tv_usec) / 1E6 + (t_end.tv_sec - t_ini.tv_sec);
}

int main(int argc, char **argv) {
    int rank, num_procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int *vres = NULL, **res = NULL;
    int *local_vres = NULL;
    int *sendcounts_rows = NULL, *displacements_rows = NULL;
    int *rows_to_compute = NULL;  // array con los índices de filas asignadas a este proceso
    int *all_rows = NULL;
    struct timeval ti, tf;

    int remainder = Y_RESN % num_procs;
    int rows_per_process = Y_RESN / num_procs;

    sendcounts_rows = (int *)malloc(num_procs * sizeof(int));
    displacements_rows = (int *)malloc(num_procs * sizeof(int));

    int displacement = 0;
    for (int i = 0; i < num_procs; i++) {
        sendcounts_rows[i] = rows_per_process + (i < remainder ? 1 : 0);
        displacements_rows[i] = displacement;
        displacement += sendcounts_rows[i];
    }

    // Solo el proceso 0 genera la lista de filas (enteros)
    if (rank == 0) {
        all_rows = (int *)malloc(Y_RESN * sizeof(int));
        for (int i = 0; i < Y_RESN; i++) {
            all_rows[i] = i;
        }
    }

    int local_row_count = sendcounts_rows[rank];
    rows_to_compute = (int *)malloc(local_row_count * sizeof(int));
    local_vres = (int *)malloc(local_row_count * X_RESN * sizeof(int));

    // ScatterV para enviar qué filas calcula cada proceso
    MPI_Scatterv(all_rows, sendcounts_rows, displacements_rows, MPI_INT,
                 rows_to_compute, local_row_count, MPI_INT,
                 0, MPI_COMM_WORLD);

    if (rank == 0) gettimeofday(&ti, NULL);

    for (int r = 0; r < local_row_count; r++) {
        int i = rows_to_compute[r];
        for (int j = 0; j < X_RESN; j++) {
            Compl z = {0.0, 0.0};
            Compl c = {
                X_MIN + j * (X_MAX - X_MIN) / X_RESN,
                Y_MAX - i * (Y_MAX - Y_MIN) / Y_RESN
            };
            int k = 0;
            float lengthsq, temp;

            do {
                temp = z.real * z.real - z.imag * z.imag + c.real;
                z.imag = 2.0 * z.real * z.imag + c.imag;
                z.real = temp;
                lengthsq = z.real * z.real + z.imag * z.imag;
                k++;
            } while (lengthsq < 4.0 && k < maxIterations);

            local_vres[r * X_RESN + j] = (k >= maxIterations) ? 0 : k;
        }
    }

    if (rank == 0) {
        vres = (int *)malloc(Y_RESN * X_RESN * sizeof(int));
        res = (int **)malloc(Y_RESN * sizeof(int *));
        for (int i = 0; i < Y_RESN; i++) {
            res[i] = vres + i * X_RESN;
        }
    }

    // Para recolectar el buffer de resultados
    int *recvcounts = (int *)malloc(num_procs * sizeof(int));
    int *displs = (int *)malloc(num_procs * sizeof(int));
    for (int i = 0; i < num_procs; i++) {
        recvcounts[i] = sendcounts_rows[i] * X_RESN;
        displs[i] = displacements_rows[i] * X_RESN;
    }

    MPI_Gatherv(local_vres, local_row_count * X_RESN, MPI_INT,
                vres, recvcounts, displs, MPI_INT,
                0, MPI_COMM_WORLD);

    if (rank == 0) {
        gettimeofday(&tf, NULL);
        fprintf(stderr, "(PERF) Time (seconds) = %lf\n", get_seconds(ti, tf));
        if (DEBUG) {
            for (int i = 0; i < Y_RESN; i++) {
                for (int j = 0; j < X_RESN; j++) {
                    printf("%3d ", res[i][j]);
                }
                printf("\n");
            }
        }
        free(vres);
        free(res);
        free(all_rows);
    }

    free(local_vres);
    free(rows_to_compute);
    free(sendcounts_rows);
    free(displacements_rows);
    free(recvcounts);
    free(displs);

    MPI_Finalize();
    return 0;
}
