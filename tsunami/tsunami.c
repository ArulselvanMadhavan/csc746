#include "file_ops.h"
#include "math.h"
#include "mpi.h"
#include "omp-test.h"
#include "omp.h"
#include "stdio.h"
#include "stdlib.h"
#include "timer.h"

#define SQ(x) ((x) * (x))
#define SWAP_PTR(xnew, xold, xtmp) (xtmp = xnew, xnew = xold, xold = xtmp)
const int ny = 200;
const int nx = 500;
const int plusYDims[] = {ny + 1, nx};
const int plusXDims[] = {ny, nx + 1};
const int nhalo = 1;
const int gdims[] = {ny + 2 * nhalo, nx + 2 * nhalo};
const double g = 9.8;
const double sigma = 0.95;
const int ntimes = 2000;
const int nburst = 100;

double **malloc2D(int jmax, int imax);

double calculateMass(double **restrict H) {
  double result = 0.0;
#pragma omp for
  for (int j = 1; j <= ny; j++) {
    for (int i = 1; i <= nx; i++) {
      result += H[j][i];
    }
  }
  return result;
}

int main3(int argc, char *argv[]) {
  test();
  return 0;
}

int main(int arc, char *argv[]) {
  const int dims[2] = {gdims[0], gdims[1]};
  const int N = gdims[0] * gdims[1];
  const int X = ny * (nx + 1);
  const int Y = (ny + 1) * nx;
  double *restrict H = malloc(sizeof(double) * N);
  double *restrict U = malloc(sizeof(double) * N);
  double *restrict V = malloc(sizeof(double) * N);
  double *restrict Hnew = malloc(sizeof(double) * N);
  double *restrict Unew = malloc(sizeof(double) * N);
  double *restrict Vnew = malloc(sizeof(double) * N);

  double *restrict Hx = malloc(sizeof(double) * X);
  double *restrict Ux = malloc(sizeof(double) * X);
  double *restrict Vx = malloc(sizeof(double) * X);
  double *restrict Hy = malloc(sizeof(double) * Y);
  double *restrict Uy = malloc(sizeof(double) * Y);
  double *restrict Vy = malloc(sizeof(double) * Y);
  double time = 0.0;
  int graph_num = 0;

  float xwin_gap = 0.0;
  float ybot_gap = 0;
  float ytop_gap = 0;
  float xwinmin = 0.0 - xwin_gap;
  float xwinmax = (float)nx + xwin_gap;
  float ywinmin = 0.0 - ybot_gap;
  float ywinmax = (float)ny + ytop_gap;
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
  /* init_io(rank, nprocs); */
  init_graphics_output();

  /* int N = gdims[0] * gdims[1]; */
  /* int X = ny * (nx + 1); */
  /* int Y = (ny + 1) * nx; */
  const double prefix = ((10.0 - 2.0) / (double)((nx + 1) / 2));
  double deltaT = 1.0e30;
  const double deltaX = 1.0, deltaY = 1.0;
  double origTM = 0;
  struct timespec starttime;
  cpu_timer_start(&starttime);
  // clang-format off
#pragma omp target data map(from:dims) map(alloc: H[:N], U[:N], V[:N], Hx[:X],Ux[:X],Vx[:X],Hy[:Y],Uy[:Y],Vy[:Y], Hnew[:N],Unew[:N],Vnew[:N]) map(to: H[:N], deltaT)
  // clang-format on
  {
    // clang-format off
#pragma omp target teams distribute parallel for collapse(2) map(always, to:dims[:2])
    // clang-format on
    for (int j = 0; j < dims[0]; j++) {
      int row = j * dims[1];
      for (int i = 0; i < dims[1]; i++) {
        int pos = row + i;
        H[pos] = 2.0;
        U[pos] = 0.0;
        V[pos] = 0.0;
      }
    }

    // clang-format off
#pragma omp target teams distribute parallel for collapse(2) map(dims[:2])
    // clang-format on
    for (int j = 0; j < dims[0]; j++) {
      int row = j * dims[1];
      for (int i = 0; i < dims[1] / 2; i++) {
        int pos = row + i;
        H[pos] = 10.0 - prefix * (double)(i);
      }
    }
    // clang-format off
/* #pragma omp target teams distribute parallel for collapse(2) map(always,from:deltaT) */
    // clang-format on
    for (int j = 1; j < ny; j++) {
      int row = j * nx;
      for (int i = 1; i < nx; i++) {
        int pos = row + i;
        double wavespeed = sqrt(g * H[pos]);
        double xspeed = (fabs(U[pos]) + wavespeed) / deltaX;
        double yspeed = (fabs(V[pos]) + wavespeed) / deltaY;
        double my_deltaT = sigma / (xspeed + yspeed);
        if (my_deltaT < deltaT)
          deltaT = my_deltaT;
      }
    }

    printf("Iteration:%d\t%d\tH:%f\tdt:%f\tMass:%f\tny:nx:%d:%d\n", dims[0],
           dims[1], H[0], deltaT, origTM, ny, nx);
    for (int n = 0; n < ntimes;) {
      for (int ib = 0; ib < nburst; ib++) {
        // clang-format off
#pragma omp target teams distribute parallel for
        // clang-format on
        for (int j = 1; j <= ny; j++) {
          int row = j * dims[1];
          H[row] = H[row + 1];
          U[row] = -U[row + 1];
          V[row] = V[row + 1];
          int last = row + nx + 1;
          H[last] = H[last - 1];
          U[last] = -U[last - 1];
          V[last] = V[last - 1];
        }

        for (int i = 0; i <= nx + 1; i++) {
          int row1 = 1 * dims[1];
          int rowN = (ny + 1) * dims[1];
          int rowNy = ny * dims[1];
          H[i] = H[row1 + i];
          U[i] = U[row1 + i];
          V[i] = -V[row1 + i];
          H[rowN + i] = H[rowNy + i];
          U[rowN + i] = U[rowNy + i];
          V[rowN + i] = -V[rowNy + i];
        }

        // clang-format off
#pragma omp target teams distribute parallel for collapse(2) map(always,from:deltaT)
        // clang-format on
        for (int j = 1; j < ny; j++) {
          int row = j * dims[1];
          for (int i = 1; i < nx; i++) {
            int pos = row + i;
            double wavespeed = sqrt(g * H[pos]);
            double xspeed = (fabs(U[pos]) + wavespeed) / deltaX;
            double yspeed = (fabs(V[pos]) + wavespeed) / deltaY;
            double my_deltaT = sigma / (xspeed + yspeed);
            if (my_deltaT < deltaT)
              deltaT = my_deltaT;
          }
        }

/* #pragma omp target map(always, from : H[:N]) */
#pragma omp target teams distribute parallel for collapse(2)
        for (int j = 0; j < ny; j++) {
          int row = j * (nx + 1);
          for (int i = 0; i <= nx; i++) {
            int pos = row + i;
            int j1i1 = (j + 1) * dims[1] + (i + 1);
            int j1i = (j + 1) * dims[1] + i;
            // density calculation
            Hx[pos] = 0.5 * (H[j1i1] + H[j1i]) -
                      deltaT / (2.0 * deltaX) * (U[j1i1] - U[j1i]);
            // momentum x calculation
            Ux[pos] = 0.5 * (U[j1i1] + U[j1i]) -
                      deltaT / (2.0 * deltaX) *
                          ((SQ(U[j1i1]) / H[j1i1] + 0.5 * g * SQ(H[j1i1])) -
                           (SQ(U[j1i]) / H[j1i] + 0.5 * g * SQ(H[j1i])));
            // momentum y calculation
            Vx[pos] =
                0.5 * (V[j1i1] + V[j1i]) - deltaT / (2.0 * deltaX) *
                                               ((U[j1i1] * V[j1i1] / H[j1i1]) -
                                                (U[j1i] * V[j1i] / H[j1i]));
          }
        }
/* #pragma omp target map(always, from : Hx[:X], Ux[:X], Vx[:X]) */

/*         printf("H:%f\tU:%f\tV:%f\n", Hx[0], Ux[240], Vx[240]); */
#pragma omp target teams distribute parallel for collapse(2)
        for (int j = 0; j <= ny; j++) {
          int row = j * nx;
          for (int i = 0; i < nx; i++) {
            int pos = row + i;
            int j1i1 = (j + 1) * dims[1] + (i + 1);
            int ji1 = j * dims[1] + (i + 1);
            // density calculation
            Hy[pos] = 0.5 * (H[j1i1] + H[ji1]) -
                      deltaT / (2.0 * deltaY) * (V[j1i1] - V[ji1]);
            // momentum x calculation
            Uy[pos] =
                0.5 * (U[j1i1] + U[ji1]) - deltaT / (2.0 * deltaY) *
                                               ((V[j1i1] * U[j1i1] / H[j1i1]) -
                                                (V[ji1] * U[ji1] / H[ji1]));
            // momentum y calculation
            Vy[pos] = 0.5 * (V[j1i1] + V[ji1]) -
                      deltaT / (2.0 * deltaY) *
                          ((SQ(V[j1i1]) / H[j1i1] + 0.5 * g * SQ(H[j1i1])) -
                           (SQ(V[ji1]) / H[ji1] + 0.5 * g * SQ(H[ji1])));
          }
        }

/* #pragma omp target map(always, from : Hy[:Y], Uy[:Y], Vy[:Y]) */

/*         printf("H:%f\tU:%f\tV:%f\n", Hy[0], Uy[240], Vy[240]); */
#pragma omp target teams distribute parallel for collapse(2)
        for (int j = 1; j <= ny; j++) {
          int row = j * dims[1];
          for (int i = 1; i <= nx; i++) {
            int pos = row + i;
            int Xj1i = (j - 1) * (nx + 1) + i;
            int Xj1i1 = (j - 1) * (nx + 1) + (i - 1);
            int Yji1 = j * nx + (i - 1);
            int Yj1i1 = j * nx + (i - 1);
            // density calculation
            Hnew[pos] = H[pos] - (deltaT / deltaX) * (Ux[Xj1i] - Ux[Xj1i1]) -
                        (deltaT / deltaY) * (Vy[Yji1] - Vy[Yj1i1]);
            // momentum x calculation
            Unew[pos] =
                U[pos] -
                (deltaT / deltaX) *
                    ((SQ(Ux[Xj1i]) / Hx[Xj1i] + 0.5 * g * SQ(Hx[Xj1i])) -
                     (SQ(Ux[Xj1i1]) / Hx[Xj1i1] + 0.5 * g * SQ(Hx[Xj1i1]))) -
                (deltaT / deltaY) * ((Vy[Yji1] * Uy[Yji1] / Hy[Yji1]) -
                                     (Vy[Yj1i1] * Uy[Yj1i1] / Hy[Yj1i1]));
            // momentum y calculation
            Vnew[pos] =
                V[pos] -
                (deltaT / deltaX) * ((Ux[Xj1i] * Vx[Xj1i] / Hx[Xj1i]) -
                                     (Ux[Xj1i1] * Vx[Xj1i1] / Hx[Xj1i1])) -
                (deltaT / deltaY) *
                    ((SQ(Vy[Yji1]) / Hy[Yji1] + 0.5 * g * SQ(Hy[Yji1])) -
                     (SQ(Vy[Yj1i1]) / Hy[Yj1i1] + 0.5 * g * SQ(Hy[Yj1i1])));
          }
        }

#pragma omp target teams distribute parallel for collapse(2)
        for (int j = 1; j <= ny; j++) {
          int row = j * dims[1];
          for (int i = 1; i <= nx; i++) {
            int pos = row + i;
            H[pos] = Hnew[pos];
            U[pos] = Unew[pos];
            V[pos] = Vnew[pos];
          }
        }
        /* printf("Iteration:%d\tdt:%f\n", n, deltaT); */
      }
      n += nburst;
      time += deltaT;
      set_data(H);
      printf("dt:%f\n", deltaT);
      /* if (rank == 0) { */
      write_to_file(graph_num, n, time);
      /* } */
      graph_num++;
    }

#pragma omp target map(always, from : H[:N])
    printf("%d\t%d\tH:%f\tdt:%f\tny:%d\n", dims[0], dims[1], H[0], deltaT, ny);
  }

  printf("Finished Target\n");
  // calculate mass
  origTM = 0;
  for (int j = 1; j <= ny; j++) {
    int row = j * dims[1];
    for (int i = 1; i <= nx; i++) {
      int pos = row + i;
      origTM += H[pos];
    }
  }
  double totalTime = cpu_timer_stop(starttime);
  printf("TotalTime:%f\tMass:%f\tdt:%f\n", totalTime, origTM, deltaT);
}

void initArrays(double **restrict H, double **restrict U, double **restrict V) {
#pragma omp for
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

int main2(int argc, char *argv[]) {
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

  /* double **restrict temp; */

  initArrays(H, U, V);

  float xwin_gap = 0.0;
  float ybot_gap = 0;
  float ytop_gap = 0;
  float xwinmin = 0.0 - xwin_gap;
  float xwinmax = (float)nx + xwin_gap;
  float ywinmin = 0.0 - ybot_gap;
  float ywinmax = (float)ny + ytop_gap;
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
  init_io(rank, nprocs);
  init_graphics_output();

  int graph_num = 0;
  /* write_to_file(graph_num, 0, 0.0); */
  /* double origTM = calculateMass(H); */

  double deltaT = 1.0e30;
  double deltaX = 1.0;
  double deltaY = 1.0;

#pragma omp for
  for (int j = 1; j < ny; j++) {
    for (int i = 1; i < nx; i++) {
      double wavespeed = sqrt(g * H[j][i]);
      double xspeed = (fabs(U[j][i]) + wavespeed) / deltaX;
      double yspeed = (fabs(V[j][i]) + wavespeed) / deltaY;
      double my_deltaT = sigma / (xspeed + yspeed);
      if (my_deltaT < deltaT) {
        deltaT = my_deltaT;
      }
    }
  }

  double time = 0.0;
  /* printf("Iteration:%5.5d, Time:%f, Timestep:%f Total mass:%f\n", 0,
   * time,
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
#pragma omp for
        for (int j = 1; j <= ny; j++) {
          H[j][0] = H[j][1];
          U[j][0] = -U[j][1]; /* why neg? */
          V[j][0] = V[j][1];
          H[j][nx + 1] = H[j][nx];
          U[j][nx + 1] = -U[j][nx];
          V[j][nx + 1] = V[j][nx];
        }

#pragma omp for
        for (int i = 0; i <= nx + 1; i++) {
          H[0][i] = H[1][i];
          U[0][i] = U[0][i];
          V[0][i] = -V[0][i];
          H[ny + 1][i] = H[ny][i];
          U[ny + 1][i] = U[ny][i];
          V[ny + 1][i] = -V[ny][i];
        }

        deltaT = 1.0e30;
#pragma omp for /* collapse */
        for (int j = 1; j < ny; j++) {
          for (int i = 1; i < nx; i++) {
            double wavespeed = sqrt(g * H[j][i]);
            double xspeed = (fabs(U[j][i]) + wavespeed) / deltaX;
            double yspeed = (fabs(V[j][i]) + wavespeed) / deltaY;
            double my_deltaT = sigma / (xspeed + yspeed);
            if (my_deltaT < deltaT) {
              deltaT = my_deltaT;
            }
          }
        }

// first pass
// x direction
#pragma omp for
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
#pragma omp for
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
#pragma omp for
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

/* SWAP_PTR(H, Hnew, temp); */
/* SWAP_PTR(U, Unew, temp); */
/* SWAP_PTR(V, Vnew, temp); */
#pragma omp for
        for (int j = 1; j <= ny; j++) {
          for (int i = 1; i <= nx; i++) {
            H[j][i] = Hnew[j][i];
            U[j][i] = Unew[j][i];
            V[j][i] = Vnew[j][i];
          }
        }
      }
    } // burst loop
    n += nburst;
    time += deltaT;
    // print iteration info
    set_data((double *)H);

    if (rank == 0) {
      double totalMass = calculateMass(H);
      printf("Iteration:%5.5d, Time:%f, Timestep:%f Mass:%f\n", n, time, deltaT,
             totalMass);
      write_to_file(graph_num, n, time);
    }
    /* printf("%d\t%p\n", rank, (double *)H); */
    /* divide_and_write(rank, gdims, graph_num, n, time); */
    /* parallel_write(graph_num, n, time); */
    graph_num++;
  }
  if (rank == 0) {
    double totaltime = cpu_timer_stop(starttime);
    printf("Rank:%d\tFlow finished in %lf seconds\n", rank, totaltime);
  }
  finalize_io();

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

  return MPI_Finalize();
}
