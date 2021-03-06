#ifndef _HDF_FILE_OPS_H
#define _HDF_FILE_OPS_H
#include "hdf5.h"
#include "mpi.h"
void hdf5_file_init(int ng, int ndims, int ny_global, int nx_global, int ny,
                    int nx, int ny_off, int nx_off, MPI_Comm mpi_hdf5_comm,
                    hid_t *memspace, hid_t *filespace);
void hdf5_file_finalize(hid_t *memspace, hid_t *filespace);

void write_hdf5_file(const char *filename, double **data, hid_t memspace,
                     hid_t filespace, MPI_Comm mpi_hdf5_comm);
void read_hdf5_file(const char *filename, double **data, hid_t memspace, hid_t filespace, MPI_Comm mpi_hdf5_comm);
#endif
