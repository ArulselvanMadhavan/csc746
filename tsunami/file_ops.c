#include "file_ops.h"
#include "hdf5_file_ops.h"
#include "memory.h"
#include "mpi.h"
#include "stddef.h"
#include "stdio.h"
#include "stdlib.h"
#include <sys/stat.h>

double *data_double = NULL;
static double xconversion = 0.0;
static double yconversion = 0.0;
const static int Ncolors = 256;
static int iteration = 0;
static int graphics_mysize = 0;
static float graphics_xmin = 0.0, graphics_xmax = 0.0, graphics_ymin = 0.0,
             graphics_ymax = 0.0;
char *graphics_directory = "graphics_output";
static double *x_double = NULL, *y_double = NULL, *dx_double = NULL,
              *dy_double = NULL;
const double scaleMax = 25.0, scaleMin = 0.0;
static const int WINSIZE = 800;

static int io_procs = 0;
static int data_per_proc = 0;
static double *data_recv;
static hid_t memspace = H5S_NULL, filespace = H5S_NULL;
static double **restrict file_buf = NULL;

double **malloc2D(int jmax, int imax);

void init_io(int rank, int nprocs) {
  io_procs = nprocs - 1;
  data_per_proc = graphics_mysize / io_procs;
  data_recv = (double *)malloc(sizeof(double) * data_per_proc);
  file_buf = malloc2D(data_per_proc, 5);
  int ng = 0;
  int ndims = 2;
  int ny_global = graphics_mysize;
  int nx_global = 5;
  int ny = data_per_proc;
  int nx = 5;

  if (rank > 0) {
    int ny_offset = (rank - 1) * data_per_proc;
    int nx_offset = 0;
    hdf5_file_init(ng, ndims, ny_global, nx_global, ny, nx, ny_offset,
                   nx_offset, MPI_COMM_WORLD, &memspace, &filespace);
  } else {
    hdf5_file_init(ng, ndims, ny_global, nx_global, 0, 0, 0, 0, MPI_COMM_WORLD,
                   &memspace, &filespace);
  }
}

void finalize_io() { hdf5_file_finalize(&memspace, &filespace); }

void set_data(double *data_in) { data_double = data_in; }
void set_graphics_mysize(int graphics_mysize_in) {
  graphics_mysize = graphics_mysize_in;
}
void set_graphics_window(float graphics_xmin_in, float graphics_xmax_in,
                         float graphics_ymin_in, float graphics_ymax_in) {
  graphics_xmin = graphics_xmin_in;
  graphics_xmax = graphics_xmax_in;
  graphics_ymin = graphics_ymin_in;
  graphics_ymax = graphics_ymax_in;
}

void set_graphics_cell_coordinates(double *x_in, double *dx_in, double *y_in,
                                   double *dy_in) {
  x_double = x_in;
  dx_double = dx_in;
  y_double = y_in;
  dy_double = dy_in;
}

void init_graphics_output() {
  /* int width = (WINSIZE / (graphics_ymax - graphics_ymin)) * */
  /* (graphics_xmax - graphics_xmin); */
  xconversion = (double)WINSIZE / (graphics_xmax - graphics_xmin);
  yconversion = (double)WINSIZE / (graphics_ymax - graphics_ymin);

  struct stat stat_descriptor;
  if (stat(graphics_directory, &stat_descriptor) == -1) {
    mkdir(graphics_directory, 0777);
  }
}

/* double** restrict data_recv = (double **)malloc2D(); */
MPI_Status status;
/* static int Ncolors = 256; */

const double step = Ncolors / (scaleMax - scaleMin);

int calc_color(double data) {
  int color = 0;
  color = (int)(data - scaleMin) * step;
  color = Ncolors - color;
  if (color < 0) {
    color = 0;
  }
  if (color >= Ncolors)
    color = Ncolors - 1;
  return color;
}

void divide_and_write(int rank, const int gdims[2], int graph_num, int ncycle,
                      double simTime) {

  MPI_Request req;
  char filename[30];
  if (rank == 0) {
    for (int i = 0; i < io_procs; i++) {
      int j_size = gdims[0];
      double *data_start = (data_double + j_size) + (i * data_per_proc);
      printf("Sending Iteration:%d\tStart:%p\n", graph_num, (void *)data_start);
      MPI_Isend(data_start, data_per_proc, MPI_DOUBLE, i + 1, graph_num,
                MPI_COMM_WORLD, &req);
    }
  } else {
    MPI_Recv(data_recv, data_per_proc, MPI_DOUBLE, 0, graph_num, MPI_COMM_WORLD,
             &status);
    printf("Received Rank:%d\n", rank);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  sprintf(filename, "graph%d.hdf5", graph_num);
  int xloc, xwid, yloc, ywid;

  if (rank > 0) {
    int start = ((rank - 1) * data_per_proc);
    int end = start + data_per_proc;
    printf("Rank:%d\tstart:%d\tend:%d\n", rank, start, end);
    for (int row = 0; row < data_per_proc; row++) {
      int i = start + row;
      xloc = (int)((x_double[i] - graphics_xmin) * xconversion);
      xwid = (int)((x_double[i] + dx_double[i] - graphics_xmin) * xconversion -
                   xloc);
      yloc =
          (int)((graphics_ymax - (y_double[i] + dy_double[i])) * yconversion);
      ywid = (int)((graphics_ymax - y_double[i]) * yconversion);
      ywid -= yloc;

      file_buf[row][0] = xloc;
      file_buf[row][1] = yloc;
      file_buf[row][2] = xwid;
      file_buf[row][3] = ywid;
      file_buf[row][4] = calc_color(data_recv[row]);
    }
  }
  printf("Rank:%d\tWriting data\n", rank);
  write_hdf5_file(filename, file_buf, memspace, filespace, MPI_COMM_WORLD);
  /* } */
  MPI_Barrier(MPI_COMM_WORLD);
}

void write_to_file(int graph_num, int ncycle, double simTime) {
  int i;
  int color;
  char filename[50], filename2[50];
  sprintf(filename, "%s/graph%05d.data", graphics_directory, graph_num);
  sprintf(filename2, "%s/outline%05d.lin", graphics_directory, graph_num);

  FILE *fp = fopen(filename, "w");
  FILE *fp2 = fopen(filename2, "w");
  if (fp && fp2) {
    fprintf(fp, "%d,%lf\n", ncycle, simTime);

    /* double step = Ncolors / (scaleMax - scaleMin); */
    int xloc, xwid, yloc, ywid;
    int xloc1, xloc2, yloc1, yloc2;
    /* printf("DISPLAY " */
    /*        "DATA\tstep:%f\tg_xmin:%f\tx_conv:%f\tg_xmax:%f\tg_ymin:%f\tg_ymax:%"
     */
    /*        "f\ts_max:%f\ts_min:%f\n", */
    /*        step, graphics_xmin, xconversion, graphics_xmax, graphics_ymin, */
    /*        graphics_ymax, scaleMax, scaleMin); */

    for (i = 0; i < graphics_mysize; i++) {
      color = calc_color(data_double[i]);
      xloc = (int)((x_double[i] - graphics_xmin) * xconversion);
      xwid = (int)((x_double[i] + dx_double[i] - graphics_xmin) * xconversion -
                   xloc);
      yloc =
          (int)((graphics_ymax - (y_double[i] + dy_double[i])) * yconversion);
      ywid = (int)((graphics_ymax - y_double[i]) * yconversion);
      ywid -= yloc;
      // fprintf(fp,"%d,%d,%d,%d,%f\n",xloc,yloc,xwid,ywid,data[i]);
      fprintf(fp, "%d,%d,%d,%d,%d\n", xloc, yloc, xwid, ywid, color);

      xloc1 = (int)((x_double[i] - graphics_xmin) * xconversion);
      xloc2 = (int)((x_double[i] + dx_double[i] - graphics_xmin) * xconversion);
      yloc1 = (int)((graphics_ymax - y_double[i]) * yconversion);
      yloc2 =
          (int)((graphics_ymax - (y_double[i] + dy_double[i])) * yconversion);
      fprintf(fp2, "%d,%d,%d,%d\n", xloc1, yloc2, xloc2, yloc2);
      fprintf(fp2, "%d,%d,%d,%d\n", xloc1, yloc1, xloc2, yloc1);
      fprintf(fp2, "%d,%d,%d,%d\n", xloc1, yloc1, xloc1, yloc2);
      fprintf(fp2, "%d,%d,%d,%d\n", xloc2, yloc1, xloc2, yloc2);
    }
    fclose(fp);
    fclose(fp2);
    iteration++;
  } else {
    if (fp == NULL) {
      printf("Could not open %s in DisplayStateToFile\n", filename);
    } else {
      printf("Could not open %s in DisplayStateToFile\n", filename2);
    }
  }
}
