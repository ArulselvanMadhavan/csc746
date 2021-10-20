//
// (C) 2021, E. Wes Bethel
// sobel_gpu.cpp
// usage:
//      sobel_gpu [no args, all is hard coded]
//

#include <chrono>
#include <iostream>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <vector>

// see https://en.wikipedia.org/wiki/Sobel_operator

// easy-to-find and change variables for the input.
// specify the name of a file containing data to be read in as bytes, along with
// dimensions [columns, rows]

// this is the original laughing zebra image
// static char input_fname[] = "../data/zebra-gray-int8";
// static int data_dims[2] = {3556, 2573};
// char output_fname[] = "../data/processed-raw-int8-gpu.dat";

// this one is a 4x augmentation of the laughing zebra
static char input_fname[] = "../data/zebra-gray-int8-4x";
static int data_dims[2] = {7112, 5146};
char output_fname[] = "../data/processed-raw-int8-4x-gpu.dat";

// see
// https://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
// macro to check for cuda errors. basic idea: wrap this macro around every cuda
// call
#define gpuErrchk(ans)                                                         \
  { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line,
                      bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
            line);
    if (abort)
      exit(code);
  }
}

__device__ float sobel_filtered_pixel(float *in, int x, int y, int rows,
                                      int cols, float *gx, float *gy) {
  int gi = 0;
  float gytemp = 0.0;
  float gxtemp = 0.0;
  for (int ky = -1; ky < 2; ky++) {
    for (int kx = -1; kx < 2; kx++) {
      int cxIdx = x + kx;
      int cyIdx = y + ky;
      int cyPos = cyIdx * cols;
      int cxy = cyPos + cxIdx;
      float sxy = in[cxy];
      gxtemp += gx[gi] * sxy;
      gytemp += gy[gi] * sxy;
      gi++;
    }
  }
  return sqrt((gxtemp * gxtemp) + (gytemp * gytemp));
}

__global__ void sobel_kernel_gpu(
    float *s, // source image pixels
    float *d, // dst image pixels
    int n,    // size of image cols*rows,
    int rows, int cols, float *gx,
    float *gy) // gx and gy are stencil weights for the sobel filter
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = gridDim.x * blockDim.x;
  for (int xy = index; xy < n; xy += stride) {
    int y = xy / cols;
    int x = xy - y * cols;
    int outIdx = y * cols + x;
    if (y == 0 || y == rows - 1 || x == 0 || x == cols - 1) {
      d[outIdx] = s[outIdx];
    } else {
      d[outIdx] = sobel_filtered_pixel(s, x, y, rows, cols, gx, gy);
    }
  }
}

int main(int ac, char *av[]) {
  // input, output file names hard coded at top of file

  // load the input file
  off_t nvalues = data_dims[0] * data_dims[1];
  unsigned char *in_data_bytes =
      (unsigned char *)malloc(sizeof(unsigned char) * nvalues);

  FILE *f = fopen(input_fname, "r");
  if (fread((void *)in_data_bytes, sizeof(unsigned char), nvalues, f) !=
      nvalues * sizeof(unsigned char)) {
    printf("Error reading input file. \n");
    fclose(f);
    return 1;
  } else
    printf(" Read data from the file %s \n", input_fname);
  fclose(f);

#define ONE_OVER_255 0.003921568627451

  // now convert input from byte, in range 0..255, to float, in range 0..1
  float *in_data_floats;
  gpuErrchk(cudaMallocManaged(&in_data_floats, sizeof(float) * nvalues));

  for (off_t i = 0; i < nvalues; i++)
    in_data_floats[i] = (float)in_data_bytes[i] * ONE_OVER_255;

  // now, create a buffer for output
  float *out_data_floats;
  gpuErrchk(cudaMallocManaged(&out_data_floats, sizeof(float) * nvalues));
  for (int i = 0; i < nvalues; i++)
    out_data_floats[i] = 1.0; // assign "white" to all output values for debug

  // define sobel filter weights, copy to a device accessible buffer
  float Gx[9] = {1.0, 0.0, -1.0, 2.0, 0.0, -2.0, 1.0, 0.0, -1.0};
  float Gy[9] = {1.0, 2.0, 1.0, 0.0, 0.0, 0.0, -1.0, -2.0, -1.0};
  float *device_gx, *device_gy;
  gpuErrchk(cudaMallocManaged(&device_gx, sizeof(float) * sizeof(Gx)));
  gpuErrchk(cudaMallocManaged(&device_gy, sizeof(float) * sizeof(Gy)));

  for (int i = 0; i < 9; i++) // copy from Gx/Gy to device_gx/device_gy
  {
    device_gx[i] = Gx[i];
    device_gy[i] = Gy[i];
  }

  // now, induce memory movement to the GPU of the data in unified memory
  // buffers

  int deviceID = 0; // assume GPU#0, always. OK assumption for this program
  cudaMemPrefetchAsync((void *)in_data_floats, nvalues * sizeof(float),
                       deviceID);
  cudaMemPrefetchAsync((void *)out_data_floats, nvalues * sizeof(float),
                       deviceID);
  cudaMemPrefetchAsync((void *)device_gx, sizeof(Gx) * sizeof(float), deviceID);
  cudaMemPrefetchAsync((void *)device_gy, sizeof(Gy) * sizeof(float), deviceID);

  // set up to run the kernel
  int nBlocks = 0, nThreadsPerBlock = 256;
  if(const char* nBlocksChar = std::getenv("NUM_BLOCKS")){
        std::cout << "Your NUM_BLOCKS is: " << nBlocksChar << '\n';
	nBlocks = std::atoi(nBlocksChar);
  } else {
      std::cout << "NUM_BLOCKS not set. Defaulting to 1" << '\n';
      nBlocks = 1;
  }
    
  // ADD CODE HERE: insert your code here to set a different number of thread
  // blocks or # of threads per block

  printf(" GPU configuration: %d blocks, %d threads per block \n", nBlocks,
         nThreadsPerBlock);

  // invoke the kernel on the device
  sobel_kernel_gpu<<<nBlocks, nThreadsPerBlock>>>(
      in_data_floats, out_data_floats, nvalues, data_dims[1], data_dims[0],
      device_gx, device_gy);

  // wait for it to finish, check errors
  gpuErrchk(cudaDeviceSynchronize());

  // write output after converting from floats in range 0..1 to bytes in range
  // 0..255
  unsigned char *out_data_bytes =
      in_data_bytes; // just reuse the buffer from before
  for (off_t i = 0; i < nvalues; i++)
    out_data_bytes[i] = (unsigned char)(out_data_floats[i] * 255.0);

  f = fopen(output_fname, "w");

  if (fwrite((void *)out_data_bytes, sizeof(unsigned char), nvalues, f) !=
      nvalues * sizeof(unsigned char)) {
    printf("Error writing output file. \n");
    fclose(f);
    return 1;
  } else
    printf(" Wrote the output file %s \n", output_fname);
  fclose(f);
}

// eof
