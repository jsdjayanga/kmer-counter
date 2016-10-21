#include <stdio.h>
#include "GPUHandler.h"

__global__ void testKernel(){
  printf("Hello from testKernel\n");
}

int64_t processKMers(char* input, char* result, int64_t size, int64_t lineLength) {
  testKernel<<<1,2>>>();
  cudaDeviceSynchronize();
  return 0;
}


