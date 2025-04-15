/*The Mandelbrot set is a fractal that is defined as the set of points c
in the complex plane for which the sequence z_{n+1} = z_n^2 + c
with z_0 = 0 does not tend to infinity.*/

/*This code computes an image of the Mandelbrot set.*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>

#define DEBUG 1

#define          X_RESN  1024  /* x resolution */
#define          Y_RESN  1024  /* y resolution */

/* Boundaries of the mandelbrot set */
#define           X_MIN  -2.0
#define           X_MAX   2.0
#define           Y_MIN  -2.0
#define           Y_MAX   2.0

/* More iterations -> more detailed image & higher computational cost */
#define   maxIterations  1000

typedef struct complextype{
  float real, imag;
} Compl;

static inline double get_seconds(struct timeval t_ini, struct timeval t_end){
  return (t_end.tv_usec - t_ini.tv_usec) / 1E6 + (t_end.tv_sec - t_ini.tv_sec);
}

int main (int argc, char *argv[]){

  /* Mandelbrot variables */
  int i, j, k;
  Compl   z, c;
  float   lengthsq, temp;
  int *vres, *res[Y_RESN];

  /* Timestamp variables */
  struct timeval  ti, tf;

  /* Allocate result matrix of Y_RESN x X_RESN */
  vres = (int *) malloc(Y_RESN * X_RESN * sizeof(int));
  if (!vres)
  {
    fprintf(stderr, "Error allocating memory\n");
    return 1;
  }
  for (i=0; i<Y_RESN; i++)
    res[i] = vres + i*X_RESN;

  /* Start measuring time */
  gettimeofday(&ti, NULL);

  /* Calculate and draw points */
  for(i=0; i < Y_RESN; i++)
  {
    for(j=0; j < X_RESN; j++)
    {
      z.real = z.imag = 0.0;
      c.real = X_MIN + j * (X_MAX - X_MIN)/X_RESN;
      c.imag = Y_MAX - i * (Y_MAX - Y_MIN)/Y_RESN;
      k = 0;

      do
      {    /* iterate for pixel color */
        temp = z.real*z.real - z.imag*z.imag + c.real;
        z.imag = 2.0*z.real*z.imag + c.imag;
        z.real = temp;
        lengthsq = z.real*z.real+z.imag*z.imag;
        k++;
      } while (lengthsq < 4.0 && k < maxIterations);

      if (k >= maxIterations) res[i][j] = 0;
      else res[i][j] = k;
    }
  }

  /* End measuring time */
  gettimeofday(&tf, NULL);
  fprintf (stderr, "(PERF) Time (seconds) = %lf\n", get_seconds(ti,tf));

  /* Print result out */
  if( DEBUG ) {
    for(i=0;i<Y_RESN;i++) {
      for(j=0;j<X_RESN;j++)
              printf("%3d ", res[i][j]);
      printf("\n");
    }
  }

  free(vres);

  return 0;
}


//CHATGPT


// #include <stdio.h>
// #include <stdlib.h>
// #include <sys/time.h>
// #include <mpi.h>

// #define          X_RESN  1024
// #define          Y_RESN  1024

// #define           X_MIN  -2.0
// #define           X_MAX   2.0
// #define           Y_MIN  -2.0
// #define           Y_MAX   2.0

// #define   maxIterations  1000

// typedef struct complextype {
//     float real, imag;
// } Compl;

// static inline double get_seconds(struct timeval t_ini, struct timeval t_end) {
//     return (t_end.tv_usec - t_ini.tv_usec) / 1E6 + (t_end.tv_sec - t_ini.tv_sec);
// }

// int main(int argc, char *argv[]) {
//     int i, j, k;
//     Compl z, c;
//     float lengthsq, temp;
//     int *local_res, *global_res;
//     struct timeval t_start, t_compute_end, t_end;
//     int rank, size;
//     int rows_per_proc;

//     MPI_Init(&argc, &argv);                      // Inicializar MPI
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);        // ID de proceso
//     MPI_Comm_size(MPI_COMM_WORLD, &size);        // Total de procesos

//     if (Y_RESN % size != 0) {
//         if (rank == 0)
//             fprintf(stderr, "Y_RESN debe ser divisible por el número de procesos\n");
//         MPI_Finalize();
//         return 1;
//     }

//     rows_per_proc = Y_RESN / size;

//     local_res = (int*) malloc(rows_per_proc * X_RESN * sizeof(int));
//     if (rank == 0) {
//         global_res = (int*) malloc(Y_RESN * X_RESN * sizeof(int));
//     }

//     MPI_Barrier(MPI_COMM_WORLD);     // Sincronizamos antes de empezar
//     gettimeofday(&t_start, NULL);

//     // Cálculo local de Mandelbrot
//     for (i = 0; i < rows_per_proc; i++) {
//         int global_i = rank * rows_per_proc + i;
//         for (j = 0; j < X_RESN; j++) {
//             z.real = z.imag = 0.0;
//             c.real = X_MIN + j * (X_MAX - X_MIN) / X_RESN;
//             c.imag = Y_MAX - global_i * (Y_MAX - Y_MIN) / Y_RESN;
//             k = 0;

//             do {
//                 temp = z.real*z.real - z.imag*z.imag + c.real;
//                 z.imag = 2.0 * z.real * z.imag + c.imag;
//                 z.real = temp;
//                 lengthsq = z.real*z.real + z.imag*z.imag;
//                 k++;
//             } while (lengthsq < 4.0 && k < maxIterations);

//             local_res[i * X_RESN + j] = (k >= maxIterations) ? 0 : k;
//         }
//     }

//     gettimeofday(&t_compute_end, NULL);

//     // Recogida de resultados
//     MPI_Gather(local_res, rows_per_proc * X_RESN, MPI_INT,
//                global_res, rows_per_proc * X_RESN, MPI_INT,
//                0, MPI_COMM_WORLD);

//     gettimeofday(&t_end, NULL);

//     // Tiempos
//     double t_compute = get_seconds(t_start, t_compute_end);
//     double t_total = get_seconds(t_start, t_end);
//     double t_comm = t_total - t_compute;

//     printf("(PERF) Rank %d - Compute: %.6lf s, Comm: %.6lf s, Total: %.6lf s\n",
//             rank, t_compute, t_comm, t_total);

//     // Solo el proceso 0 imprime la imagen completa (opcional)
//     #ifdef DEBUG
//     if (rank == 0) {
//         for (i = 0; i < Y_RESN; i++) {
//             for (j = 0; j < X_RESN; j++)
//                 printf("%3d ", global_res[i * X_RESN + j]);
//             printf("\n");
//         }
//     }
//     #endif

//     free(local_res);
//     if (rank == 0) free(global_res);

//     MPI_Finalize();
//     return 0;
// }
