//
// (C) 2021, E. Wes Bethel
// mpi_2dmesh.cpp

// usage:
//      mpi_2dmesh [args as follows]
//
// command line arguments:
//    -g 1|2|3  : domain decomp method, where 1=row-slab decomp, 2=column-slab
//    decomp,
//          3=tile decomp (OPTIONAL, default is -g 1, row-slab decomp)
//    -i filename : name of datafile containing raw unsigned bytes
//          (OPTIONAL, default input filename in mpi_2dmesh.hpp)
//    -x XXX : specify the number of columns in the mesh, the width (REQUIRED)
//    -y YYY : specify the number of rows in the mesh, the height (REQUIRED)
//    -o filename : where output results will be written (OPTIONAL, default
//    output
//          filename in mpi_2dmesh.hpp)
//    -a 1|2 : the action to perform: 1 means perform per-tile processing then
//    write
//       output showing results, 2 means generate an output file with cells
//       labeled as to which rank they belong to (depends on -g 1|2|3 setting)
//       (OPTIONAL, default: -a 1, perform the actual processing)
//    -v : a flag that will trigger printing out the 2D vector array of Tile2D
//    (debug)
//       (OPTIONAL, default value is no verbose debug output)
//
// Assumptions:
//
// Grid decompositions:
//       When creating tile-based decomps, will compute the number of tiles in
//       each dimension as sqrt(nranks); please use a number of ranks that is
//       the square of an integer, e.g., 4, 9, 16, 25, etc.

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <unistd.h>
#include <vector>

#include "mpi_2dmesh.hpp" // for AppState and Tile2D class

#define DEBUG_TRACE 0
    int nhalo = 1;
int parseArgs(int ac, char *av[], AppState *as) {
  int rstat = 0;
  int c;
  // Arul: No output filename?
  while ((c = getopt(ac, av, "va:g:x:y:i:")) != -1) {
    switch (c) {
    case 'a': {
      int action = std::atoi(optarg == NULL ? "-1" : optarg);
      if (action != MESH_PROCESSING && action != MESH_LABELING_ONLY) {
        printf("Error parsing command line: -a %s is an undefined action \n",
               optarg);
        rstat = 1;
      } else
        as->action = action;
      break;
    }

    case 'g': {
      int decomp = std::atoi(optarg == NULL ? "-1" : optarg);
      if (decomp != ROW_DECOMP && decomp != COLUMN_DECOMP &&
          decomp != TILE_DECOMP) {
        printf(
            "Error parsing command line: -g %s is an undefined decomposition\n",
            optarg);
        rstat = 1;
      } else
        as->decomp = decomp;
      break;
    }
    case 'x': {
      int xval = std::atoi(optarg == NULL ? "-1" : optarg);
      if (xval < 0) {
        printf(" Error parsing command line: %d is not a valid mesh width \n",
               xval);
        rstat = 1;
      } else
        as->global_mesh_size[0] = xval;
      break;
    }
    case 'y': {
      int yval = std::atoi(optarg == NULL ? "-1" : optarg);
      if (yval < 0) {
        printf(" Error parsing command line: %d is not a valid mesh height \n",
               yval);
        rstat = 1;
      } else
        as->global_mesh_size[1] = yval;
      break;
    }
    case 'i': {
      strcpy(as->input_filename, optarg);
      break;
    }
    case 'o': {
      strcpy(as->output_filename, optarg);
      break;
    }
    case 'v': {
      as->debug = 1;
      break;
    }
    } // end switch

  } // end while

  return rstat;
}

//
// computeMeshDecomposition:
// input: AppState *as - nranks, mesh size, etc.
// computes: tileArray, a 2D vector of Tile2D objects
//
// assumptions:
// - For tile decompositions, will use sqrt(nranks) tiles per axis
//
void computeMeshDecomposition(AppState *as, vector<vector<Tile2D>> *tileArray) {
  int xtiles, ytiles;
  int ntiles;

  if (as->decomp == ROW_DECOMP) {
    // in a row decomposition, each tile width is the same as the mesh width
    // the mesh is decomposed along height
    xtiles = 1;

    //  set up the y coords of the tile boundaries
    ytiles = as->nranks;
    int ylocs[ytiles + 1];
    int ysize = as->global_mesh_size[1] / as->nranks; // size of each tile in y

    int yval = 0;
    for (int i = 0; i < ytiles; i++, yval += ysize) {
      ylocs[i] = yval;
    }
    ylocs[ytiles] = as->global_mesh_size[1];

    // then, create tiles along the y axis
    for (int i = 0; i < ytiles; i++) {
      vector<Tile2D> tiles;
      int width = as->global_mesh_size[0];
      int height = ylocs[i + 1] - ylocs[i];
      Tile2D t = Tile2D(0, ylocs[i], width, height, i, 0, i, 1, ytiles);
      tiles.push_back(t);
      tileArray->push_back(tiles);
    }
  } else if (as->decomp == COLUMN_DECOMP) {
    // in a columne decomposition, each tile height is the same as the mesh
    // height the mesh is decomposed along width
    ytiles = 1;

    // set up the x coords of the tile boundaries
    xtiles = as->nranks;
    int xlocs[xtiles + 1];
    int xsize = as->global_mesh_size[0] / as->nranks; // size of each tile in x

    int xval = 0;
    for (int i = 0; i < xtiles; i++, xval += xsize) {
      xlocs[i] = xval;
    }
    xlocs[xtiles] = as->global_mesh_size[0];

    // then, create tiles along the x axis
    vector<Tile2D> tile_row;
    for (int i = 0; i < xtiles; i++) {
      int width = xlocs[i + 1] - xlocs[i];
      int height = as->global_mesh_size[1];
      Tile2D t = Tile2D(xlocs[i], 0, width, height, i, i, 0, xtiles, 1);
      tile_row.push_back(t);
    }
    tileArray->push_back(tile_row);
  } else // assume as->decom == TILE_DECOMP
  {
    // to keep things simple, we  will assume sqrt(nranks) tiles in each of x
    // and y axes. if sqrt(nranks) is not an even integer, then this approach
    // will result in some ranks without tiles/work to do

    double root = sqrt(as->nranks);
    int nranks_per_axis = (int)root;

    xtiles = ytiles = nranks_per_axis;

    // set up x coords for tile boundaries
    int xlocs[xtiles + 1];
    int xsize =
        as->global_mesh_size[0] / nranks_per_axis; // size of each tile in x

    int xval = 0;
    for (int i = 0; i < xtiles; i++, xval += xsize) {
      xlocs[i] = xval;
    }
    xlocs[xtiles] = as->global_mesh_size[0];

    // set up y coords for tile boundaries
    int ylocs[ytiles + 1];
    int ysize =
        as->global_mesh_size[1] / nranks_per_axis; // size of each tile in y

    int yval = 0;
    for (int i = 0; i < ytiles; i++, yval += ysize) {
      ylocs[i] = yval;
    }
    ylocs[ytiles] = as->global_mesh_size[1];

    // now, build 2D array of tiles
    int rank = 0;
    for (int j = 0; j < ytiles; j++) { // fix me
      vector<Tile2D> tile_row;
      for (int i = 0; i < xtiles; i++) {
        int width, height;
        width = xlocs[i + 1] - xlocs[i];
        height = ylocs[j + 1] - ylocs[j];
        Tile2D t = Tile2D(xlocs[i], ylocs[j], width, height, rank++, i, j,
                          xtiles, ytiles);
        tile_row.push_back(t);
      }
      tileArray->push_back(tile_row);
    }
  }
}

void write_output_labels(AppState as, vector<vector<Tile2D>> tileArray) {
  // create a buffer of ints, we will set buf[i,j] to be a value
  // in the range 0..nranks-1 to reflect which rank is owner of
  // a particular buf[i,j] location.

  size_t xsize = as.global_mesh_size[0];
  size_t ysize = as.global_mesh_size[1];
  size_t gsize = xsize * ysize;
  printf("\n\nWriting out mesh labels to a file \n");

  int *meshRankLabels = (int *)malloc(gsize * sizeof(int));

  for (off_t i = 0; i < xsize * ysize; i++)
    meshRankLabels[i] = -1;

  for (int row = 0; row < tileArray.size(); row++) {
    printf(" Row %d of the tileArray is of length %zu \n", row,
           tileArray[row].size());

    for (int col = 0; col < tileArray[row].size(); col++) {
      Tile2D t = tileArray[row][col];

      t.print(row, col); // row, col

      // for each tile, paint the meshRankLabels array with the tile's rank
      // at the mesh locations owned by the tile

      // base offset into the output buffer for this tile
      size_t bufOffset = t.yloc * xsize + t.xloc;

      for (int jj = 0; jj < t.height; jj++, bufOffset += xsize) {
        for (int ii = 0; ii < t.width; ii++)
          meshRankLabels[bufOffset + ii] = t.tileRank;
      }

      xsize = as.global_mesh_size[0]; // breakpoint anchor
    }                                 // end loop over tileArray columns
    ysize = as.global_mesh_size[1];   // breakpoint anchor
  }                                   // end loop over tileArray rows

  FILE *f = fopen(as.output_filename, "w");
  if (f == NULL) {
    perror(" Error opening output file ");
  }
  fwrite((void *)meshRankLabels, sizeof(int), gsize, f);
  fclose(f);
  free(meshRankLabels);
} // end writing an output buffer

void printTileArray(vector<vector<Tile2D>> &tileArray) {
  printf("---- Contents of the tileArray, which is of length %zu \n",
         tileArray.size());

  for (int row = 0; row < tileArray.size(); row++) {
    printf(" Row %d of the tileArray is of length %zu \n", row,
           tileArray[row].size());
    for (int col = 0; col < tileArray[row].size(); col++) {
      Tile2D t = tileArray[row][col];
      t.print(row, col);
    }
  }
}

float byteNormalize(unsigned char i) {
#define ONE_OVER_255 0.003921568627451
  return ((float)i * ONE_OVER_255);
}

void loadInputFile(AppState *as) {
  // open the input file
  FILE *f = fopen(as->input_filename, "r");
  if (f == NULL) {
    perror(" Problem opening input file ");
    return;
  }

  // create the landing buffer for the input data
  size_t nitems = as->global_mesh_size[0] * as->global_mesh_size[1];
  vector<unsigned char> buf(nitems);

  // read the input into the buffer
  size_t bytes_read =
      fread((void *)buf.data(), sizeof(unsigned char), nitems, f);
  fclose(f);

  // create space for the float-converted buffer
  as->input_data_floats.resize(nitems);

  // now convert from byte, in range 0..255, to float, in range 0..1
  transform(buf.begin(), buf.end(), as->input_data_floats.begin(),
            byteNormalize);
}

unsigned char floatNormalize(float t) {
  // assume t in range 0..1
  return ((unsigned char)(t * 255.0));
}

void writeOutputFile(AppState &as) {
  // open the input file
  FILE *f = fopen(as.output_filename, "w");
  if (f == NULL) {
    perror(" Problem opening output file ");
    return;
  }

#if 0
   // this code writes out the output as floats rather than converting to ubyte
   size_t nitems = as.global_mesh_size[0] * as.global_mesh_size[1];
   fwrite((const void *)as.output_data_floats.data(), sizeof(float), nitems, f);
#endif

#if 1
  // create the landing buffer for the input data
  // size_t nitems = as.global_mesh_size[1] * as.global_mesh_size[0];
  size_t nitems = as.output_data_floats.size();
  vector<unsigned char> buf(nitems);

  // now convert from byte, in range 0..255, to float, in range 0..1
  transform(as.output_data_floats.begin(), as.output_data_floats.end(),
            buf.begin(), floatNormalize);
  // float *out = as.output_data_floats.data();
  // int Gw = as.global_mesh_size[0] + 2 * nhalo;
  // int Gh = as.global_mesh_size[1] + 2 * nhalo;

  // int out_off = (nhalo * Gw) + nhalo;
  // for (int j = 0; j < as.global_mesh_size[1]; j++) {
  //   int b_off = j * as.global_mesh_size[0];
  //   int o_off = out_off + (j * Gw);
  //   // printf("b_off:%d\to_off:%d\n", b_off, o_off);
  //   for (int i = 0; i < as.global_mesh_size[0]; i++) {
  //     int bPos = b_off + i;
  //     int oPos = o_off + i;

  //     buf[bPos] = floatNormalize(out[oPos]);
  //   }
  // }
  // write out the byte buffer
  fwrite((const void *)buf.data(), sizeof(unsigned char), nitems, f);
#endif

  fclose(f);
}

void gsendStridedBuffer(float *srcBuf, int srcW, int srcH, int sendW, int sendH,
                        int fromRank, int toRank) {
  int msgTag = 0;
  int bufcount = sendW * sendH;
  float *buffer = (float *)malloc(bufcount * sizeof(float));
  int d_offset = 0;
  int s_offset = (nhalo * srcW) + nhalo;
  for (int j = 0; j < sendH; j++, d_offset += sendW, s_offset += srcW) {
    memcpy((void *)(buffer + d_offset), (void *)(srcBuf + s_offset),
           (sizeof(float) * sendW));
  }
  MPI_Send(buffer, bufcount, MPI_FLOAT, toRank, msgTag, MPI_COMM_WORLD);
}

void sendStridedBuffer(float *srcBuf, int srcWidth, int srcHeight,
                       int srcOffsetColumn, int srcOffsetRow, int sendWidth,
                       int sendHeight, int fromRank, int toRank) {
  int msgTag = 0;
  MPI_Request request;
  int srcPos = srcOffsetRow * srcWidth + srcOffsetColumn;
  int count = sendWidth * sendHeight;
  float *startPos = srcBuf + srcPos;
  float *buffer = (float *)malloc(count * sizeof(float));
  for (int j = 0; j < sendHeight; j++) {
    int bPos = j * sendWidth;
    int rPos = j * srcWidth;
    memcpy((void *)(buffer + bPos), (void *)(startPos + rPos),
           sizeof(float) * sendWidth);
  }
  MPI_Send(buffer, count, MPI_FLOAT, toRank, msgTag, MPI_COMM_WORLD);
  // MPI_Request_free(&request);
}

void recvStridedBuffer(float *dstBuf, int dstWidth, int dstHeight,
                       int dstOffsetColumn, int dstOffsetRow, int expectedWidth,
                       int expectedHeight, int fromRank, int toRank) {

  int msgTag = 0;
  MPI_Status status[4];
  MPI_Request request[4];

  int ecount = expectedHeight * expectedWidth;
  float *buffer = (float *)malloc(ecount * sizeof(float));
  MPI_Recv(buffer, ecount, MPI_FLOAT, fromRank, msgTag, MPI_COMM_WORLD,
           &status[0]);
  int s_off = 0;
  int Dw = dstWidth + 2 * nhalo;
  int d_off = ((dstOffsetRow + nhalo) * Dw) + (dstOffsetColumn + nhalo);
  // printf("dr:%d\tdc:%d\tstart:%d\n", dstOffsetRow + nhalo,
  // dstOffsetColumn + nhalo, d_off);
  for (int j = 0; j < expectedHeight;
       j++, s_off += expectedWidth, d_off += Dw) {
    // printf("Soff:%d\tdoff:%d\teW:%d\n", s_off, d_off, expectedWidth);
    memcpy((void *)(dstBuf + d_off), (void *)(buffer + s_off),
           sizeof(float) * expectedWidth);
  }
}

float sobel_filtered_pixel(float *in, int x, int y, int dims[], float *gx,
                           float *gy) {
  int gi = 0;
  float sum = 0.0;
  float gytemp = 0.0;
  float gxtemp = 0.0;
  for (int ky = -1; ky < 2; ky++) {
    for (int kx = -1; kx < 2; kx++) {
      int cxIdx = x + kx;
      int cyIdx = y + ky;
      int cyPos = cyIdx * dims[0];
      int cxy = cyPos + cxIdx;
      float sxy = in[cxy];
      gxtemp += gx[gi] * sxy;
      gytemp += gy[gi] * sxy;
      gi++;
    }
  }
  return sqrt((gxtemp * gxtemp) + (gytemp * gytemp));
}

// Rows * columns
float Gx[] = {1.0, 0.0, -1.0, 2.0, 0.0, -2.0, 1.0, 0.0, -1.0};
float Gy[] = {1.0, 2.0, 1.0, 0.0, 0.0, 0.0, -1.0, -2.0, -1.0};

void sobelAllTiles(int myrank, vector<vector<Tile2D>> &tileArray) {
  for (int row = 0; row < tileArray.size(); row++) {
    for (int col = 0; col < tileArray[row].size(); col++) {
      Tile2D *t = &(tileArray[row][col]);

      if (t->tileRank == myrank) {
        // std::copy(t->inputBuffer.begin(), t->inputBuffer.end(),
        // t->outputBuffer.begin());
        float *in = t->inputBuffer.data();
        float *out = t->outputBuffer.data();
        int dims[2] = {t->width, t->height};
        int ghdims[2] = {t->width + 2 * nhalo, t->height + 2 * nhalo};
        // for (int y = 1; y < dims[1] - 1; y++) {
        //   for (int x = 1; x < dims[0] - 1; x++) {

        //     int outIdx = y * dims[0] + x;
        //     // out[outIdx] = in[outIdx];
        //     out[outIdx] = sobel_filtered_pixel(in, x, y, dims, Gx, Gy);
        //   }
        // }
        for (int y = 0; y < ghdims[1]; y++) {
          int off = y * ghdims[0];
          for (int x = 0; x < ghdims[0]; x++) {
            int outIdx = off + x;
            out[outIdx] = in[outIdx];
          }
        }
      }
    }
  }
}

void scatterAllTiles(int myrank, vector<vector<Tile2D>> &tileArray, float *s,
                     int global_width, int global_height) {

#if DEBUG_TRACE
  printf(" Rank %d is entering scatterAllTiles \n", myrank);
#endif
  for (int row = 0; row < tileArray.size(); row++) {
    for (int col = 0; col < tileArray[row].size(); col++) {
      Tile2D *t = &(tileArray[row][col]);
      int tileW = t->width + 2 * nhalo;
      int tileH = t->height + 2 * nhalo;
      int tileSize = tileW * tileH;
      if (myrank != 0 && t->tileRank == myrank) {
        int fromRank = 0;

        // receive a tile's buffer
        t->inputBuffer.resize(tileSize, 0);
        t->outputBuffer.resize(tileSize, 0);
#if DEBUG_TRACE
        printf("scatterAllTiles() receive side:: t->tileRank=%d, myrank=%d, "
               "t->inputBuffer->size()=%zu, t->outputBuffersize()=%zu \n",
               t->tileRank, myrank, t->inputBuffer.size(),
               t->outputBuffer.size());
#endif
        recvStridedBuffer(
            t->inputBuffer.data(), t->width, t->height, 0,
            0, // offset into the tile buffer: we want the whole thing
            t->width, t->height, // how much data coming from this tile
            fromRank, myrank);
      } else if (myrank == 0) {
        if (t->tileRank != 0) {
#if DEBUG_TRACE
          printf("scatterAllTiles() send side: t->tileRank=%d, myrank=%d, "
                 "t->inputBuffer->size()=%zu \n",
                 t->tileRank, myrank, t->inputBuffer.size());
#endif
          sendStridedBuffer(s, // ptr to the buffer to send
                            global_width,
                            global_height,       // size of the src buffer
                            t->xloc, t->yloc,    // offset into the send buffer
                            t->width, t->height, // size of the buffer to send,
                            myrank, t->tileRank);
        } else // rather then have rank 0 send to rank 0, just do a strided copy
               // into a tile's input buffer
        {
          t->inputBuffer.resize(tileSize, 0);
          t->outputBuffer.resize(tileSize, 0);

          off_t s_offset = 0, d_offset = (nhalo * tileW) + nhalo;
          float *d = t->inputBuffer.data();

          for (int j = 0; j < t->height;
               j++, s_offset += global_width, d_offset += tileW) {
            memcpy((void *)(d + d_offset), (void *)(s + s_offset),
                   sizeof(float) * t->width);
          }
        }
      }
    }
  } // loop over 2D array of tiles

#if DEBUG_TRACE
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 1) {
    printf("\n\n ----- rank=%d, inside scatterAllTiles debug printing of the "
           "tile array \n",
           myrank);
    printTileArray(tileArray);
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void gatherAllTiles(int myrank, vector<vector<Tile2D>> &tileArray, float *d,
                    int global_width, int global_height) {
  for (int row = 0; row < tileArray.size(); row++) {
    for (int col = 0; col < tileArray[row].size(); col++) {
      Tile2D *t = &(tileArray[row][col]);
      int tileW = t->width + 2 * nhalo;
      int tileH = t->height + 2 * nhalo;
      int tileSize = tileW * tileH;
#if DEBUG_TRACE
      printf("gatherAllTiles(): t->tileRank=%d, myrank=%d, "
             "t->outputBuffer->size()=%zu \n",
             t->tileRank, myrank, t->outputBuffer.size());
#endif

      if (myrank != 0 && t->tileRank == myrank) {
        gsendStridedBuffer(t->outputBuffer.data(), tileW, tileH, t->width,
                           t->height, t->tileRank, 0);
      } else if (myrank == 0) {
        if (t->tileRank != 0) {
          // receive a tile's buffer and copy back into the output buffer d
          // printf("Collecting output(buffered output):%d(%dx%d)\n", tileSize,
          // tileW, tileH);
          recvStridedBuffer(d, global_width, global_height, t->xloc,
                            t->yloc, // offset of this tile
                            t->width,
                            t->height, // how much data coming from this tile
                            t->tileRank, myrank);
        } else // copy from a tile owned by rank 0 back into the main buffer
        {
          float *s = t->outputBuffer.data();
          int Gw = global_width + 2 * nhalo;
          off_t s_offset = (nhalo * tileW) + nhalo,
                d_offset = ((t->yloc + nhalo) * Gw) + (t->xloc + nhalo);

          for (int j = 0; j < t->height;
               j++, s_offset += tileW, d_offset += Gw) {
            memcpy((void *)(d + d_offset), (void *)(s + s_offset),
                   sizeof(float) * t->width);
          }
        }
      }
    }
  } // loop over 2D array of tiles
}

void ghost_update(int myRank, vector<vector<Tile2D>> &tileArray) {
  // Lookup scatter impl
  // Iterate over 2d tile vector
  for (int row = 0; row < tileArray.size(); row++) {
    for (int col = 0; col < tileArray[row].size(); col++) {
      Tile2D *t = &(tileArray[row][col]);
      if (t->tileRank == myRank) {
        printf("%d\tUpdating ghost\n", myRank);
      }
    }
  }
}

int main(int ac, char *av[]) {

  AppState as;
  vector<vector<Tile2D>> tileArray;

  std::chrono::time_point<std::chrono::high_resolution_clock> start_time,
      end_time;
  std::chrono::duration<double> elapsed_scatter_time, elapsed_ghost_time,
      elapsed_sobel_time, elapsed_gather_time;

  MPI_Init(&ac, &av);

  int myrank, nranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  as.myrank = myrank;
  as.nranks = nranks;

  if (parseArgs(ac, av, &as) != 0) {
    MPI_Finalize();
    return 1;
  }

  char hostname[256];
  gethostname(hostname, sizeof(hostname));

  printf("\nHello world, I'm rank %d of %d total ranks running on <%s>\n",
         as.myrank, as.nranks, hostname);
  MPI_Barrier(MPI_COMM_WORLD);

#if DEBUG_TRACE
  if (as.myrank == 0)
    printf("\n\n ----- All ranks will computeMeshDecomposition \n");
#endif

  computeMeshDecomposition(&as, &tileArray);

  if (as.myrank == 0 && as.debug == 1) // print out the AppState and tileArray
  {
    as.print();
    printTileArray(tileArray);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (as.action == MESH_LABELING_ONLY) {
    if (as.myrank == 0) {
      printf("\n\n Rank 0 is writing out mesh labels to a file \n");
      write_output_labels(as, tileArray);
    }
  } else {
    // === rank 0 loads the input file
    if (as.myrank == 0) {
#if DEBUG_TRACE
      printf("\n\n Rank 0 is loading input \n");
#endif
      loadInputFile(&as);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // ----------- scatter phase of processing

    // start the timer
    start_time = std::chrono::high_resolution_clock::now();

    scatterAllTiles(as.myrank, tileArray, as.input_data_floats.data(),
                    as.global_mesh_size[0], as.global_mesh_size[1]);

    // end the timer
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = std::chrono::high_resolution_clock::now();
    elapsed_scatter_time = end_time - start_time;

    // ----------- the actual processing
    MPI_Barrier(MPI_COMM_WORLD);

    start_time = std::chrono::high_resolution_clock::now();
    ghost_update(myrank, tileArray);
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = std::chrono::high_resolution_clock::now();
    elapsed_ghost_time = end_time - start_time;

    // start the timer
    start_time = std::chrono::high_resolution_clock::now();

    sobelAllTiles(as.myrank, tileArray);

    // end the timer
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = std::chrono::high_resolution_clock::now();
    elapsed_sobel_time = end_time - start_time;

    // ----------- gather processing
    int Gw = as.global_mesh_size[0] + 2 * nhalo;
    int Gh = as.global_mesh_size[1] + 2 * nhalo;
    // create output buffer space on rank 0
    if (as.myrank == 0) {
      as.output_data_floats.resize(Gw * Gh);
      // initialize to a known value outside the range of expected values
      std::fill(as.output_data_floats.begin(), as.output_data_floats.end(),
                -1.0);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // start the timer
    start_time = std::chrono::high_resolution_clock::now();

    gatherAllTiles(as.myrank, tileArray, as.output_data_floats.data(),
                   as.global_mesh_size[0], as.global_mesh_size[1]);

    // end the timer
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = std::chrono::high_resolution_clock::now();
    elapsed_gather_time = end_time - start_time;

    // === write output

    if (as.myrank == 0) // only rank 0 writes data to disk
    {
      writeOutputFile(as);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (as.myrank == 0) {
    printf("\n\nTiming results from rank 0: \n");
    printf("\tScatter time:\t%6.4f (ms) \n",
           (elapsed_scatter_time * 1000.0).count());
    printf("\tGhost time:\t%6.4f (ms) \n",
           (elapsed_ghost_time * 1000.0).count());
    printf("\tSobel time:\t%6.4f (ms) \n",
           (elapsed_sobel_time * 1000.0).count());
    printf("\tGather time:\t%6.4f (ms) \n",
           (elapsed_gather_time * 1000.0).count());
  }

  MPI_Finalize();
  return 0;
}
// EOF
