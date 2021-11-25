#include "file_ops.h"
#include "stddef.h"
#include "stdio.h"
#include <sys/stat.h>

static double *data_double = NULL;
static double xconversion = 0.0;
static double yconversion = 0.0;
static int Ncolors = 256;
static int iteration = 0;
static int graphics_mysize = 0;
static float graphics_xmin = 0.0, graphics_xmax = 0.0, graphics_ymin = 0.0,
             graphics_ymax = 0.0;
char *graphics_directory = "graphics_output";
static double *x_double = NULL, *y_double = NULL, *dx_double = NULL,
              *dy_double = NULL;
const double scaleMax = 25.0, scaleMin = 0.0;
static const int WINSIZE = 800;

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
  int width = (WINSIZE / (graphics_ymax - graphics_ymin)) *
              (graphics_xmax - graphics_xmin);
  xconversion = (double)WINSIZE / (graphics_xmax - graphics_xmin);
  yconversion = (double)WINSIZE / (graphics_ymax - graphics_ymin);

  struct stat stat_descriptor;
  if (stat(graphics_directory, &stat_descriptor) == -1) {
    mkdir(graphics_directory, 0777);
  }
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

    double step = Ncolors / (scaleMax - scaleMin);
    int xloc, xwid, yloc, ywid;
    int xloc1, xloc2, yloc1, yloc2;
    /* printf("DISPLAY " */
    /*        "DATA\tstep:%f\tg_xmin:%f\tx_conv:%f\tg_xmax:%f\tg_ymin:%f\tg_ymax:%" */
    /*        "f\ts_max:%f\ts_min:%f\n", */
    /*        step, graphics_xmin, xconversion, graphics_xmax, graphics_ymin, */
    /*        graphics_ymax, scaleMax, scaleMin); */

    for (i = 0; i < graphics_mysize; i++) {
      color = (int)(data_double[i] - scaleMin) * step;
      color = Ncolors - color;
      if (color < 0) {
        color = 0;
      }
      if (color >= Ncolors)
        color = Ncolors - 1;

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
