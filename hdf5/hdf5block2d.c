#include "hdf5.h"
#include "hdf5_file_ops.h"
#include "malloc2D.h"
#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"

void init_array(int ny, int nx, int ng, double **array) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for (int j = 0; j < ny + 2 * ng; j++) {
    for (int i = 0; i < nx + 2 * ng; i++) {
      array[j][i] = 0.0;
    }
  }

  int icount = 1;
  for (int j = ng; j < ny + ng; j++) {
    for (int i = ng; i < nx + ng; i++) {
      array[j][i] = (double)(icount + 100 * rank);
      icount++;
    }
  }
}

const int primaryRank = 0;
const int primaryColor = 0;

void assign_hdf5_ranks(int rank, int nprocs, MPI_Comm comm,
                       MPI_Comm *mpi_hdf5_comm) {
  int nfiles = 1;
  /* int ranks_per_file = nprocs - 1; */
  int color = primaryColor;
  if (rank != primaryRank) {
    color = 1;
  }
  /* color = rank % 2; */
  MPI_Comm_split(comm, color, rank, mpi_hdf5_comm);
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int rank, nprocs;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  MPI_Comm mpi_hdf5_comm = MPI_COMM_NULL;
  assign_hdf5_ranks(rank, nprocs, comm, &mpi_hdf5_comm);

  int hdf5_nprocs;
  MPI_Comm_size(mpi_hdf5_comm, &hdf5_nprocs);
  printf("Rank:%d\tHDF5 size:%d\n", rank, hdf5_nprocs);
  int ntimes = 4;
  int ng = 0;
  int ndims = 2;
  int ny_global = 4;
  int nx_global = 10;
  int ny = 1;
  int nx = 10;
  int ny_offset = rank;
  int nx_offset = 0;
  hid_t memspace = H5S_NULL, filespace = H5S_NULL;
  hdf5_file_init(ng, ndims, ny_global, nx_global, ny, nx, ny_offset, nx_offset,
                 mpi_hdf5_comm, &memspace, &filespace);
  char filename[30];
  double **data = (double **)malloc2D(ny_global, nx_global);
  MPI_Status stat;
  int rows_per_rank = ny_global / nprocs;
  for (int gnum = 0; gnum < ntimes; gnum++) {
    sprintf(filename, "example%d.hdf5", gnum);
    if (rank == 0) {
      for (int j = 0; j < ny_global; j++) {
        for (int i = 0; i < nx_global; i++) {
          data[j][i] = (gnum * 100) + (j * 10) + i;
        }
      }

      // Generate rows_per_rank
      int recv_rank = 0;
      for (int j = 1; j < ny_global; j++) {
        recv_rank = j;
        MPI_Send((double *)(data[j]), (rows_per_rank * nx_global), MPI_DOUBLE,
                 recv_rank, gnum, comm);
      }
    } else {
      double **recv_data = (double **)malloc2D(rows_per_rank, nx_global);
      MPI_Recv((double *)(recv_data[rows_per_rank - 1]),
               (rows_per_rank * nx_global), MPI_DOUBLE, 0, gnum, comm, &stat);
      write_hdf5_file(filename, recv_data, memspace, filespace, mpi_hdf5_comm);
    }
    MPI_Barrier(comm);
  }
  /* MPI_status */
  hdf5_file_finalize(&memspace, &filespace);
  MPI_Finalize();
}

/* SPMD style */
int main2(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int rank, nprocs;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  /* Split */
  MPI_Comm mpi_hdf5_comm = MPI_COMM_NULL;
  int nfiles = 1;
  float ranks_per_file = (float)nprocs / (float)nfiles;
  int color = (int)((float)rank / (float)ranks_per_file);
  MPI_Comm_split(comm, color, rank, &mpi_hdf5_comm);
  /* printf("Rank:%d/%d\tColor/RPF:%d/%f\n", rank, nprocs, color,
   * ranks_per_file); */

  int nprocs_color, rank_color;
  MPI_Comm_rank(mpi_hdf5_comm, &rank_color);
  MPI_Comm_size(mpi_hdf5_comm, &nprocs_color);
  /* printf("ColorId:%d\tRC:%d\tNC:%d\n", color, rank_color, nprocs_color); */

  MPI_Comm row_comm, col_comm;
  int row_color = 1, col_color = rank_color;
  MPI_Comm_split(mpi_hdf5_comm, row_color, rank_color, &row_comm);
  MPI_Comm_split(mpi_hdf5_comm, col_color, rank_color, &col_comm);

  int ndims = 2, ng = 2, ny = 10, nx = 10, ny_offset = 0, nx_offset = 0,
      nx_global = 0, ny_global = 0;
  int global_subsizes[] = {ny, nx};
  MPI_Exscan(&nx, &nx_offset, 1, MPI_INT, MPI_SUM, row_comm);
  MPI_Exscan(&ny, &ny_offset, 1, MPI_INT, MPI_SUM, col_comm);
  printf("rank:%d\trank_c:%d\tnx_off:%d\tny_off:%d\n", rank, rank_color,
         nx_offset, ny_offset);
  int row_size, col_size;
  MPI_Comm_size(row_comm, &row_size);
  MPI_Comm_size(col_comm, &col_size);
  printf("Row size:%d\tCol size:%d\n", row_size, col_size);

  MPI_Allreduce(&nx, &nx_global, 1, MPI_INT, MPI_SUM, row_comm);
  MPI_Allreduce(&ny, &ny_global, 1, MPI_INT, MPI_SUM, col_comm);
  printf("G_off\tRank:%d,nx:%d\tny:%d\n", rank, nx_global, ny_global);

  double **data = (double **)malloc2D(ny + 2 * ng, nx + 2 * ng);
  double **data_restore = (double **)malloc2D(ny + 2 * ng, nx + 2 * ng);

  init_array(ny, nx, ng, data);
  for (int j = 0; j < ny + 2 * ng; j++) {
    for (int i = 0; i < nx + 2 * ng; i++) {
      data_restore[j][i] = 0.0;
    }
  }

  hid_t memspace = H5S_NULL, filespace = H5S_NULL;
  hdf5_file_init(ng, ndims, ny_global, nx_global, ny, nx, ny_offset, nx_offset,
                 mpi_hdf5_comm, &memspace, &filespace);

  char filename[30];
  if (nfiles > 1) {
    sprintf(filename, "example_%02d.hdf5", color);
  } else {
    sprintf(filename, "example.hdf5");
  }

  write_hdf5_file(filename, data, memspace, filespace, mpi_hdf5_comm);
  read_hdf5_file(filename, data_restore, memspace, filespace, mpi_hdf5_comm);
  hdf5_file_finalize(&memspace, &filespace);

  if (rank == 0)
    printf("Verifying  checkpoint\n");

  int ierr = 0;
  // verification
  for (int j = 0; j < ny + 2 * ng; j++) {
    for (int i = 0; i < nx + 2 * ng; i++) {
      if (data_restore[j][i] != data[j][i]) {
        ierr++;
        printf("DEBUG -- j %d i %d restored %lf data %lf\n", j, i,
               data_restore[j][i], data[j][i]);
      }
    }
  }
  int ierr_global = 0;
  MPI_Allreduce(&ierr, &ierr_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0 && ierr_global == 0)
    printf("   Checkpoint has been verified\n");
  free(data);
  free(data_restore);

  MPI_Comm_free(&mpi_hdf5_comm);
  MPI_Comm_free(&row_comm);
  MPI_Comm_free(&col_comm);
  MPI_Finalize();
}
