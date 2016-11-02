#include <stdio.h>
#include <inttypes.h>
#include "GPUHandler.h"

__global__ void bitEncode(char* input, char* filter, int64_t lineLength) {
	uint64_t index = (blockIdx.x * blockDim.x + threadIdx.x) * lineLength;

	uint64_t readValue = 0;
	uint64_t filterValue = 0;

	uint16_t readValueLength = 0;
	uint16_t filterValueLength = 0;

	int64_t i = index;
	for (; i < index + lineLength; i++) {
		if (i - index > 0 && (i - index) % 32 == 0) {
			int64_t readValueLocation = ((((i - index) / 32) - 1)
					* sizeof(int64_t)) + sizeof(int16_t) + index;
			//int64_t filterLocation = ((((i - index) / 32) - 1) * sizeof (int64_t));
			memcpy(&input[readValueLocation], &readValue, sizeof(int64_t));
			//memcpy(&filter[filterLocation], &filterValue, sizeof (int64_t));
			readValue = 0;

			//filterValue = 0;
		}

		if (i > 0 && (i - index) % 64 == 0) {
			int64_t filterLocation =
					((((i - index) / 64) - 1) * sizeof(int64_t)) + index;
			memcpy(&filter[filterLocation], &filterValue, sizeof(int64_t));
			filterValue = 0;
		}

		switch (input[i]) {
		case 'A':
			readValue <<= 2;
			readValue |= 0;
			readValueLength += 2;

			filterValue <<= 1;
			filterValueLength += 1;

			continue;
		case 'C':
			readValue <<= 2;
			readValue |= 1;
			readValueLength += 2;

			filterValue <<= 1;
			filterValueLength += 1;

			continue;
		case 'G':
			readValue <<= 2;
			readValue |= 2;
			readValueLength += 2;

			filterValue <<= 1;
			filterValueLength += 1;

			continue;
		case 'T':
			readValue <<= 2;
			readValue |= 3;
			readValueLength += 2;

			filterValue <<= 1;
			filterValueLength += 1;

			continue;
		default:
			readValue <<= 2;
			readValue |= 3;
			readValueLength += 2;

			filterValue <<= 1;
			filterValue |= 1;
			filterValueLength += 1;
			continue;
		}
	}

	//printf("readValueLength==============:%"PRIu16"", readValueLength);
	memcpy(&input[index], &readValueLength, sizeof(uint16_t));

	if (i > 0 && (i - index) % 64 > 0) {
		uint8_t shiftingReadValue = (32 - ((i - index) % 32)) * 2;
		readValue <<= shiftingReadValue;

		int64_t readValueLocation = ((((i - index) / 32)) * sizeof(int64_t))
				+ sizeof(int16_t) + index;
		memcpy(&input[readValueLocation], &readValue, sizeof(int64_t));
		readValue = 0;
		readValueLength = 0;

		uint8_t shiftingFilterValue = 64 - ((i - index) % 64);
		filterValue <<= shiftingFilterValue;
		int64_t filterLocation = ((((i - index) / 64)) * sizeof(int64_t))
				+ index;
		memcpy(&filter[filterLocation], &filterValue, sizeof(int64_t));
		filterValue = 0;
		filterValueLength = 0;
	}

}

int64_t processKMers(const char* input, int64_t kmerLength, int64_t inputSize,
		int64_t lineLength) {
	printf("Processing k-mers klen=%"PRIu64", inSize=%"PRIu64","
	" liLen=%"PRIu64"\n", kmerLength, inputSize, lineLength);
	bool debug = true;

	char* d_input;
	char* d_output;
	char* d_filter;

	cudaMalloc((void **) &d_input, inputSize);
	cudaMalloc((void **) &d_output, inputSize);
	cudaMalloc((void **) &d_filter, inputSize / 2);

	cudaMemcpy(d_input, input, inputSize, cudaMemcpyHostToDevice);
	cudaMemset(d_output, 0, inputSize);
	cudaMemset(d_filter, 0, inputSize / 2);

	int32_t threadCount = 256;
	int32_t count = inputSize / lineLength / threadCount;
	if ((inputSize / lineLength) % threadCount > 0) {
		count++;
	}

	for (int32_t ite = 0; ite <= count; ite++) {
		bitEncode<<<1, threadCount>>>(&d_input[threadCount * lineLength * ite], d_filter, lineLength);
		cudaDeviceSynchronize();
	}

	if (debug == true) {
		char* temp = new char[inputSize];
		memset(temp, 0, inputSize);

		cudaMemcpy(temp, d_input, inputSize, cudaMemcpyDeviceToHost);

		printf(
				"%"PRIu16" : %"PRIu64", %"PRIu64", %"PRIu64", %"PRIu64", %"PRIu64"\n",
				*(uint16_t*) &temp[0], *(uint64_t*) &temp[2],
				*(uint64_t*) &temp[10], *(uint64_t*) &temp[18],
				*(uint64_t*) &temp[26], *(uint64_t*) &temp[34]);

		for (int i = 0; i < inputSize; i += lineLength) {
			uint16_t* count = (uint16_t*) &temp[i];
			printf("===============Count:%"PRIu16"  %i\n", *count, i);
//			int iterations = (*count) / 64;
//			if ((*count) % 64 > 0) {
//				iterations++;
//			}
//			for (int j = 2; j < iterations * 8; j += 8) {
//				printf("%d : %"PRIu64"\n", j, *((uint64_t*) (&temp[i + j])));
//			}
		}

		for (int i = 0; i < inputSize; i++) {
			printf("%c", temp[i]);
		}
	}

	cudaDeviceReset();

	return 0;
}

