//
// (C) 2021, E. Wes Bethel
// conv_harness_cpu.cpp:
// usage:
//      conv_harness_cpu [no args, all is hard coded]
//

#include <chrono>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <unistd.h>
#include <vector>

// easy-to-find and change variables for the input.
// specify the name of a file containing data to be read in as bytes, along with
// dimensions [columns, rows]

// this is the original laughing zebra image
// static char input_fname[] = "../data/zebra-gray-int8";
// static int data_dims[2] = {3556, 2573};
// char output_fname[] = "../data/processed-raw-int8-cpu.dat";

// this one is a 4x augmentation of the laughing zebra
static char input_fname[] = "../data/zebra-gray-int8-4x";
// static int data_dims[2] = {6, 5};
static int data_dims[2] = {7112, 5146};
static int paddingWidth = 0;
static int padded_dims[2] = {data_dims[0] + paddingWidth,
                             data_dims[1] + paddingWidth};
// static int padded_dims[2] = {data_dims[0] + 2, data_dims[1] + 2};
char output_fname[] = "../data/processed-raw-int8-4x-cpu.dat";

//
// sobel_filtered_pixel(): perform the sobel filtering at a given i,j location
//
// input: float *s - the source data
// input: int i,j - the location of the pixel in the source data where we want
// to center our sobel convolution input: int rows, cols: the dimensions of the
// input and output image buffers input: float *gx, gy:  arrays of length 9
// each, these are logically 3x3 arrays of sobel filter weights
//
// this routine computes Gx=gx*s centered at (i,j), Gy=gy*s centered at (i,j),
// and returns G = sqrt(Gx^2 + Gy^2)

// see https://en.wikipedia.org/wiki/Sobel_operator
//
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
float Gx[] = {
  1.0, 0.0, -1.0,
  2.0, 0.0, -2.0,
  1.0, 0.0, -1.0
};
float Gy[] = {
  1.0, 2.0, 1.0,
  0.0, 0.0, 0.0,
  -1.0, -2.0, -1.0
};

void do_sobel_filtering(float *in, float *out, int dims[2]) {

#pragma omp parallel default(none) shared(in, out, dims, Gx, Gy)
  {
#pragma omp for collapse(2)
    for (int y = 1; y < dims[1] - 1; y++) {
      for (int x = 1; x < dims[0] - 1; x++) {
        int outIdx = y * dims[0] + x;
        out[outIdx] = sobel_filtered_pixel(in, x, y, dims, Gx, Gy);
      }
    }
  }
}

int main(int ac, char *av[]) {

  off_t nvalues = data_dims[0] * data_dims[1];
  off_t padded_n = padded_dims[0] * padded_dims[1];
  unsigned char *in_data_bytes =
      (unsigned char *)malloc(sizeof(unsigned char) * nvalues);

  FILE *f = fopen(input_fname, "r");
  if (f == NULL) {
    printf(" Error opening the input file: %s \n", input_fname);
    return 1;
  }
  if (fread((void *)in_data_bytes, sizeof(unsigned char), nvalues, f) !=
      nvalues * sizeof(unsigned char)) {
    printf("Error reading input file. \n");
    fclose(f);
    return 1;
  } else
    printf(" Read data from the file %s \n", input_fname);
  fclose(f);

#define ONE_OVER_255 0.003921568627451

  // now convert from byte, in range 0..255, to float, in range 0..1
  float *in_data_floats = (float *)malloc(sizeof(float) * nvalues);
  for (off_t i = 0; i < nvalues; i++)
    in_data_floats[i] = (float)in_data_bytes[i] * ONE_OVER_255;

  // now, create a buffer for output
  float *out_data_floats = (float *)malloc(sizeof(float) * padded_n);
  for (int i = 0; i < padded_n; i++) {
    out_data_floats[i] = 1.0f;
  }
  // do the processing =======================
  std::chrono::time_point<std::chrono::high_resolution_clock> start_time =
      std::chrono::high_resolution_clock::now();

  do_sobel_filtering(in_data_floats, out_data_floats, data_dims);

  std::chrono::time_point<std::chrono::high_resolution_clock> end_time =
      std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = end_time - start_time;
  std::cout << " Elapsed time is : " << elapsed.count() << " " << std::endl;

  // write output after converting from floats in range 0..1 to bytes in range
  // 0..255
  unsigned char *out_data_bytes =
      (unsigned char *)malloc(sizeof(unsigned char) * padded_n);
  for (off_t i = 0; i < padded_n; i++)
    out_data_bytes[i] = (unsigned char)(out_data_floats[i] * 255.0);

  f = fopen(output_fname, "w");

  if (fwrite((void *)out_data_bytes, sizeof(unsigned char), padded_n, f) !=
      padded_n * sizeof(unsigned char)) {
    printf("Error writing output file. \n");
    fclose(f);
    return 1;
  } else
    printf(" Wrote the output file %s \n", output_fname);
  fclose(f);
}

// eof
