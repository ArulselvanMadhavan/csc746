#include "file_ops.h"
#include "math.h"
#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "timer.h"

#define SQ(x) ((x) * (x))
#define SWAP_PTR(xnew, xold, xtmp) (xtmp = xnew, xnew = xold, xold = xtmp)

const int ny = 200;
const int nx = 500;
const int nhalo = 1;
const double g = 9.8;
const double sigma = 0.95;
const int ntimes = 2000;
const int nburst = 100;
const int gdims[] = {ny + 2 * nhalo, nx + 2 * nhalo};

double **malloc2D(int jmax, int imax);

void initArrays(double **restrict H, double **restrict U, double **restrict V) {
  for (int j = 0; j < gdims[0]; j++) {
    for (int i = 0; i < gdims[1]; i++) {
      H[j][i] = 2.0;
      U[j][i] = 0.0;
      V[j][i] = 0.0;
    }

    const double prefix = ((10.0 - 2.0) / (double)((nx + 1) / 2));
    for (int i = 0; i < gdims[1] / 2; i++) {
      H[j][i] = 10.0 - prefix * (double)(i);
    }
  }
}

void calculateDeltaT(double **restrict H, double **restrict U,
                     double **restrict V, double *deltaX, double *deltaY,
                     double *deltaT) {

  for (int j = 1; j < ny; j++) {
    for (int i = 1; i < nx; i++) {
      double wavespeed = sqrt(g * H[j][i]);
      double xspeed = (fabs(U[j][i]) + wavespeed) / *deltaX;
      double yspeed = (fabs(V[j][i]) + wavespeed) / *deltaY;
      double my_deltaT = sigma / (xspeed + yspeed);
      if (my_deltaT < *deltaT) {
        *deltaT = my_deltaT;
      }
    }
  }
}

double calculateMass(double **restrict H) {
  double result = 0.0;
  for (int j = 1; j <= ny; j++) {
    for (int i = 1; i <= nx; i++) {
      result += H[j][i];
    }
  }
  return result;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int rank, nprocs;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  printf("Running Shallow Water equations with rank:%d\tnprocs:%d\n", rank,
         nprocs);
  /* Height */
  double **restrict H = malloc2D(gdims[0], gdims[1]);
  /* Mass flux = velocity * mass */
  double **restrict U = malloc2D(gdims[0], gdims[1]);
  double **restrict V = malloc2D(gdims[0], gdims[1]);

  double **restrict Hnew = malloc2D(gdims[0], gdims[1]);
  double **restrict Unew = malloc2D(gdims[0], gdims[1]);
  double **restrict Vnew = malloc2D(gdims[0], gdims[1]);

  double **restrict Hx = malloc2D(ny, nx + 1);
  double **restrict Ux = malloc2D(ny, nx + 1);
  double **restrict Vx = malloc2D(ny, nx + 1);

  double **restrict Hy = malloc2D(ny + 1, nx);
  double **restrict Uy = malloc2D(ny + 1, nx);
  double **restrict Vy = malloc2D(ny + 1, nx);

  double **restrict temp;

  initArrays(H, U, V);

  float xwinmin = 0.0 - 2.0;
  float xwinmax = (float)nx + 2.0;
  float ywinmin = 0.0 - 12.0;
  float ywinmax = (float)ny + 2.0;
  set_graphics_mysize(nx * ny);
  set_graphics_window(xwinmin, xwinmax, ywinmin, ywinmax);
  double **restrict dx = malloc2D(ny + 2, nx + 2);
  double **restrict dy = malloc2D(ny + 2, nx + 2);
  double **restrict x = malloc2D(ny + 2, nx + 2);
  double **restrict y = malloc2D(ny + 2, nx + 2);

  for (int j = 0; j <= ny + 1; j++) {
    for (int i = 0; i <= nx + 1; i++) {
      dx[j][i] = 1.0;
      dy[j][i] = 1.0;
      x[j][i] = 0.0 + (double)i * 1.0;
      y[j][i] = 0.0 + (double)j * 1.0;
    }
  }
  set_graphics_cell_coordinates((double *)x, (double *)dx, (double *)y,
                                (double *)dy);
  set_data((double *)H);
  init_graphics_output();

  int graph_num = 0;
  write_to_file(graph_num, 0, 0.0);
  double origTM = calculateMass(H);

  double deltaT = 1.0e30;
  double deltaX = 1.0;
  double deltaY = 1.0;
  calculateDeltaT(H, U, V, &deltaX, &deltaY, &deltaT);
  double time = 0.0;
  /* printf("Iteration:%5.5d, Time:%f, Timestep:%f Total mass:%f\n", 0, time,
   */
  /* deltaT, origTM); */

  struct timespec starttime;
  cpu_timer_start(&starttime);

  /* ntimes * nburst */
  for (int n = 0; n < ntimes;) {
    /* printf("Outer:%d\n", n); */
    for (int ib = 0; ib < nburst; ib++) {
      if (rank == 0) {
        /* Boundary conditions */
        for (int j = 1; j <= ny; j++) {
          H[j][0] = H[j][1];
          U[j][0] = -U[j][1]; /* why neg? */
          V[j][0] = V[j][1];
          H[j][nx + 1] = H[j][nx];
          U[j][nx + 1] = -U[j][nx];
          V[j][nx + 1] = V[j][nx];
        }

        for (int i = 0; i <= nx + 1; i++) {
          H[0][i] = H[1][i];
          U[0][i] = U[0][i];
          V[0][i] = -V[0][i];
          H[ny + 1][i] = H[ny][i];
          U[ny + 1][i] = U[ny][i];
          V[ny + 1][i] = -V[ny][i];
        }

        deltaT = 1.0e30;
        calculateDeltaT(H, U, V, &deltaX, &deltaY, &deltaT);

        // first pass
        // x direction
        for (int j = 0; j < ny; j++) {
          for (int i = 0; i <= nx; i++) {
            // density calculation
            Hx[j][i] =
                0.5 * (H[j + 1][i + 1] + H[j + 1][i]) -
                deltaT / (2.0 * deltaX) * (U[j + 1][i + 1] - U[j + 1][i]);
            // momentum x calculation
            Ux[j][i] = 0.5 * (U[j + 1][i + 1] + U[j + 1][i]) -
                       deltaT / (2.0 * deltaX) *
                           ((SQ(U[j + 1][i + 1]) / H[j + 1][i + 1] +
                             0.5 * g * SQ(H[j + 1][i + 1])) -
                            (SQ(U[j + 1][i]) / H[j + 1][i] +
                             0.5 * g * SQ(H[j + 1][i])));
            // momentum y calculation
            Vx[j][i] =
                0.5 * (V[j + 1][i + 1] + V[j + 1][i]) -
                deltaT / (2.0 * deltaX) *
                    ((U[j + 1][i + 1] * V[j + 1][i + 1] / H[j + 1][i + 1]) -
                     (U[j + 1][i] * V[j + 1][i] / H[j + 1][i]));
          }
        }

        // y direction
        for (int j = 0; j <= ny; j++) {
          for (int i = 0; i < nx; i++) {
            // density calculation
            Hy[j][i] =
                0.5 * (H[j + 1][i + 1] + H[j][i + 1]) -
                deltaT / (2.0 * deltaY) * (V[j + 1][i + 1] - V[j][i + 1]);
            // momentum x calculation
            Uy[j][i] =
                0.5 * (U[j + 1][i + 1] + U[j][i + 1]) -
                deltaT / (2.0 * deltaY) *
                    ((V[j + 1][i + 1] * U[j + 1][i + 1] / H[j + 1][i + 1]) -
                     (V[j][i + 1] * U[j][i + 1] / H[j][i + 1]));
            // momentum y calculation
            Vy[j][i] = 0.5 * (V[j + 1][i + 1] + V[j][i + 1]) -
                       deltaT / (2.0 * deltaY) *
                           ((SQ(V[j + 1][i + 1]) / H[j + 1][i + 1] +
                             0.5 * g * SQ(H[j + 1][i + 1])) -
                            (SQ(V[j][i + 1]) / H[j][i + 1] +
                             0.5 * g * SQ(H[j][i + 1])));
          }
        }

        // second pass
        for (int j = 1; j <= ny; j++) {
          for (int i = 1; i <= nx; i++) {
            // density calculation
            Hnew[j][i] = H[j][i] -
                         (deltaT / deltaX) * (Ux[j - 1][i] - Ux[j - 1][i - 1]) -
                         (deltaT / deltaY) * (Vy[j][i - 1] - Vy[j - 1][i - 1]);
            // momentum x calculation
            Unew[j][i] =
                U[j][i] -
                (deltaT / deltaX) * ((SQ(Ux[j - 1][i]) / Hx[j - 1][i] +
                                      0.5 * g * SQ(Hx[j - 1][i])) -
                                     (SQ(Ux[j - 1][i - 1]) / Hx[j - 1][i - 1] +
                                      0.5 * g * SQ(Hx[j - 1][i - 1]))) -
                (deltaT / deltaY) *
                    ((Vy[j][i - 1] * Uy[j][i - 1] / Hy[j][i - 1]) -
                     (Vy[j - 1][i - 1] * Uy[j - 1][i - 1] / Hy[j - 1][i - 1]));
            // momentum y calculation
            Vnew[j][i] =
                V[j][i] -
                (deltaT / deltaX) *
                    ((Ux[j - 1][i] * Vx[j - 1][i] / Hx[j - 1][i]) -
                     (Ux[j - 1][i - 1] * Vx[j - 1][i - 1] / Hx[j - 1][i - 1])) -
                (deltaT / deltaY) * ((SQ(Vy[j][i - 1]) / Hy[j][i - 1] +
                                      0.5 * g * SQ(Hy[j][i - 1])) -
                                     (SQ(Vy[j - 1][i - 1]) / Hy[j - 1][i - 1] +
                                      0.5 * g * SQ(Hy[j - 1][i - 1])));
          }
        }

        SWAP_PTR(H, Hnew, temp);
        SWAP_PTR(U, Unew, temp);
        SWAP_PTR(V, Vnew, temp);
      }
    } // burst loop
    n += nburst;
    time += deltaT;
    double TotalMass = calculateMass(H);
    // print iteration info
    /* printf("Iteration:%5.5d, Time:%f, Timestep:%f Total mass:%f\n", n,
     * time, */
    /* deltaT, TotalMass); */
    set_data((double *)H);
    if (rank == 0) {
      write_to_file(graph_num, n, time);
    }
    parallel_write(graph_num, n, time);
    /* MPI_Barrier(comm); */
    graph_num++;
  }

  double totaltime = cpu_timer_stop(starttime);
  printf(" Flow finished in %lf seconds\n", totaltime);

  free(H);
  free(U);
  free(V);
  free(Hnew);
  free(Unew);
  free(Vnew);
  free(Hx);
  free(Ux);
  free(Vx);
  free(Hy);
  free(Uy);
  free(Vy);
  free(x);
  free(y);
  free(dx);
  free(dy);

  MPI_Finalize();
}
