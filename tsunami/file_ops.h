#ifndef _FILE_OPS_H_
#define _FILE_OPS_H_
void set_data(double *data_in);
void init_graphics_output();
void set_graphics_cell_coordinates(double *x_in, double *dx_in, double *y_in,
                                   double *dy_in);
void set_graphics_window(float graphics_xmin_in, float graphics_xmax_in,
                         float graphics_ymin_in, float graphics_ymax_in);
void set_graphics_mysize(int graphics_mysize_in);
void write_to_file(int graph_num, int ncycle, double simTime);
void parallel_write(int graph_num, int ncycle, double simTime);
#endif
