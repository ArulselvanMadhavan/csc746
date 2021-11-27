#include "hdf5_file_ops.h"
#include "H5Dpublic.h"
#include "H5FDmpi.h"
#include "H5FDmpio.h"
#include "H5Fpublic.h"
#include "H5Ipublic.h"
#include "H5Ppublic.h"
#include "H5Spublic.h"
#include "H5Tpublic.h"
#include "H5public.h"
#include "mpi.h"

hid_t create_filespace(int ndims, int ny_global, int nx_global, int ny, int nx,
                       int ny_off, int nx_off, MPI_Comm mpi_hdf5_comm) {
  hsize_t dims[] = {ny_global, nx_global};
  hid_t filespace = H5Screate_simple(ndims, dims, NULL);

  hsize_t start[] = {ny_off, nx_off};
  hsize_t stride[] = {1, 1};
  hsize_t count[] = {ny, nx};

  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, stride, count, NULL);
  return filespace;
}

hid_t create_memspace(int ndims, int ny, int nx, int ng) {
  hsize_t dims[] = {ny + 2 * ng, nx + 2 * ng};
  hid_t ms = H5Screate_simple(ndims, dims, NULL);
  hsize_t start[] = {ng, ng};
  hsize_t stride[] = {1, 1};
  hsize_t count[] = {ny, nx};
  H5Sselect_hyperslab(ms, H5S_SELECT_SET, start, stride, count, NULL);
  return ms;
}

void hdf5_file_init(int ng, int ndims, int ny_global, int nx_global, int ny,
                    int nx, int ny_off, int nx_off, MPI_Comm mpi_hdf5_comm,
                    hid_t *memspace, hid_t *filespace) {
  *filespace = create_filespace(ndims, ny_global, nx_global, ny, nx, ny_off,
                                nx_off, mpi_hdf5_comm);
  *memspace = create_memspace(ndims, ny, nx, ng);
}

void hdf5_file_finalize(hid_t *memspace, hid_t *filespace) {
  H5Sclose(*memspace);
  *memspace = H5S_NULL;
  H5Sclose(*filespace);
  *filespace = H5S_NULL;
}

hid_t create_hdf5_file(const char *filename, MPI_Comm mpi_hdf5_comm) {
  hid_t fc_plist = H5P_DEFAULT;
  hid_t fa_plist = H5P_DEFAULT;
  fa_plist = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_coll_metadata_write(fa_plist, true);

  MPI_Info mpi_info = MPI_INFO_NULL;
  MPI_Info_create(&mpi_info);
  MPI_Info_set(mpi_info, "striping_factor", "8");
  MPI_Info_set(mpi_info, "striping_unit", "4194304");

  H5Pset_fapl_mpio(fa_plist, mpi_hdf5_comm, mpi_info);

  hid_t fid = H5Fcreate(filename, H5F_ACC_TRUNC, fc_plist, fa_plist);
  H5Pclose(fa_plist);
  MPI_Info_free(&mpi_info);
  return fid;
}

hid_t create_dataset(hid_t file_id, hid_t filespace) {
  hid_t link_creation_plist = H5P_DEFAULT;
  hid_t ds_creation_plist = H5P_DEFAULT;
  hid_t ds_access_plist = H5P_DEFAULT;
  hid_t dataset =
      H5Dcreate2(file_id, "data array", H5T_IEEE_F64LE, filespace,
                 link_creation_plist, ds_creation_plist, ds_access_plist);
  return dataset;
}

void write_hdf5_file(const char *filename, double **data, hid_t memspace,
                     hid_t filespace, MPI_Comm mpi_hdf5_comm) {
  hid_t fid = create_hdf5_file(filename, mpi_hdf5_comm);
  hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
  hid_t ds = create_dataset(fid, filespace);
  H5Dwrite(ds, H5T_IEEE_F64LE, memspace, filespace, xfer_plist, &(data[0][0]));
  H5Dclose(ds);
  H5Pclose(xfer_plist);
  H5Fclose(fid);
}

hid_t open_hdf5_file(const char *filename, MPI_Comm mpi_hdf5_comm) {
  hid_t fa_plist = H5P_DEFAULT;
  fa_plist = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_all_coll_metadata_ops(fa_plist, true);
  H5Pset_fapl_mpio(fa_plist, mpi_hdf5_comm, MPI_INFO_NULL);
  hid_t fid = H5Fopen(filename, H5F_ACC_RDONLY, fa_plist);
  H5Pclose(fa_plist);
  return fid;
}

hid_t open_hdf5_dataset(hid_t fid) {
  hid_t da_plist = H5P_DEFAULT;
  hid_t dataset = H5Dopen2(fid, "data array", da_plist);
  return dataset;
}

void read_hdf5_file(const char *filename, double **data, hid_t memspace,
                    hid_t filespace, MPI_Comm mpi_hdf5_comm) {
  hid_t fid = open_hdf5_file(filename, mpi_hdf5_comm);
  hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);

  hid_t ds = open_hdf5_dataset(fid);
  H5Dread(ds, H5T_IEEE_F64LE, memspace, filespace, H5P_DEFAULT, &(data[0][0]));
  H5Dclose(ds);
  H5Pclose(xfer_plist);
  H5Fclose(fid);
}
