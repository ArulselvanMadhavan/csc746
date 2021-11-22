#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "timer.h"

const int ny = 200;
const int nx = 500;
const int nhalo = 1;
const double g = 9.8;
const double sigma = 0.95;
const int ntimes = 2000;
const int nburst = 100;

double **malloc2D(int jmax, int imax);

int main(int argc, char *argv[]) {
  printf("Running Shallow Water equations");
  int gdims[] = {ny + 2 * nhalo, nx + 2 * nhalo};
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

  for (int j = 0; j < gdims[0]; j++) {
    for (int i = 0; i < gdims[1]; i++) {
      H[j][i] = 2.0;
      U[j][i] = 0.0;
      V[j][i] = 0.0;
    }

    const double prefix = ((10.0 - 2.0) / (double)((nx + 1) / 2));
    for (int i = 0; i <= (nx + 1) / 2; i++) {
      H[j][i] = 10.0 - prefix * (double)(i);
    }
  }

  double origTM = 0.0;
  for (int j = 1; j <= ny; j++) {
    for (int i = 1; i <= nx; i++) {
      origTM += H[j][i];
    }
  }

  double deltaT = 1.0e30;
  double deltaX = 1.0;
  double deltaY = 1.0;
  for (int j = 1; j < ny; j++) {
    for (int i = 1; i < nx; i++) {
      double wavespeed = sqrt(g * H[j][i]); /* Why sqrt */
      double xspeed = (fabs(U[j][i]) + wavespeed) / deltaX;
      double yspeed = (fabs(U[j][i]) + wavespeed) / deltaY;
      double my_deltaT = sigma / (xspeed + yspeed);
      if (my_deltaT < deltaT)
        deltaT = my_deltaT;
    }
  }
  double time = 0.0;
  printf("Iteration:%5.5d, Time:%f, Timestep:%f Total mass:%f\n", 0, time,
         deltaT, origTM);

  struct timespec starttime;
  cpu_timer_start(&starttime);

  /* ntimes * nburst */
  for (int n = 0; n < ntimes;) {
    /* printf("Outer:%d\n", n); */
    for (int ib = 0; ib < nburst; ib++) {
      /* printf("%d\t%d\n", n, ib); */
      n += nburst;
    }
  }

  free(H);
  free(U);
  free(V);
  return 0;
}
