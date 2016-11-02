#include <stdio.h>
#include <inttypes.h>
#include "GPUHandler.h"

__global__ void bitEncode(char* input, char* filter, int64_t lineLength,
		int64_t upperBound) {
	uint64_t index = (blockIdx.x * blockDim.x + threadIdx.x) * lineLength;

	if (index > upperBound - lineLength) {
		return;
	}

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

		if (i - index > 0 && (i - index) % 64 == 0) {
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

__device__ bool checkBit(uint64_t filter, uint8_t bit) {
	//cout << "++++++++++:" << (int32_t)bit << endl;
	uint64_t t = 1;
	uint64_t temp = t << (64 - bit - 1);
	if ((temp & filter) > 0) {
		return true;
	}
	return false;
}

__device__ uint64_t read64bits(char* input, int64_t index) {
	uint64_t value = 0;
	memcpy(&value, &input[index], sizeof(int64_t));
	return value;
}

__global__ void extractKMers(char* input, char* bitFilter, char*output,
		uint64_t sectionLength, int64_t kmerLength) {
	uint64_t index = (blockIdx.x * blockDim.x + threadIdx.x) * sectionLength;

}

int64_t processKMers(const char* input, int64_t kmerLength, int64_t inputSize,
		int64_t lineLength) {
	printf("Processing k-mers klen=%"PRIu64", inSize=%"PRIu64","
	" liLen=%"PRIu64"\n", kmerLength, inputSize, lineLength);
	bool debug = true;

	char* d_input;
	char* d_output;
	char* d_filter;

	uint64_t outputSize = (inputSize / lineLength) * (kmerLength / 4)
			* (lineLength - kmerLength + 1);

	cudaMalloc((void **) &d_input, inputSize);
	cudaMalloc((void **) &d_output, outputSize);
	cudaMalloc((void **) &d_filter, inputSize);

	cudaMemcpy(d_input, input, inputSize, cudaMemcpyHostToDevice);
	cudaMemset(d_output, 0, outputSize);
	cudaMemset(d_filter, 0, inputSize);

	int32_t threadCount = 256;
	int32_t count = inputSize / lineLength / threadCount;
	if ((inputSize / lineLength) % threadCount > 0) {
		count++;
	}

	for (int32_t ite = 0; ite <= count; ite++) {
		bitEncode<<<1, threadCount>>>(&d_input[threadCount * lineLength * ite],
				&d_filter[threadCount * lineLength * ite], lineLength,
				inputSize);
		cudaDeviceSynchronize();

		extractKMers<<<1, threadCount>>>(
				&d_input[threadCount * lineLength * ite],
				&d_filter[threadCount * lineLength * ite],
				&d_output[threadCount * outputSize / (inputSize / lineLength)
						* ite], outputSize / (inputSize / lineLength),
				kmerLength);
		cudaDeviceSynchronize();
	}

	printBitEncodedResult(d_input, d_filter, inputSize, lineLength);

	cudaDeviceReset();

	return 0;
}

void printBitEncodedResult(char* d_input, char* d_filter, uint64_t inputSize,
		uint64_t lineLength) {
	char* temp = new char[inputSize];
	memset(temp, 0, inputSize);

	char* tempFilter = new char[inputSize];
	memset(tempFilter, 0, inputSize);

	cudaMemcpy(temp, d_input, inputSize, cudaMemcpyDeviceToHost);
	cudaMemcpy(tempFilter, d_filter, inputSize, cudaMemcpyDeviceToHost);

	printf(
			"%"PRIu16" : %"PRIu64", %"PRIu64", %"PRIu64", %"PRIu64", %"PRIu64"\n",
			*(uint16_t*) &temp[0], *(uint64_t*) &temp[2],
			*(uint64_t*) &temp[10], *(uint64_t*) &temp[18],
			*(uint64_t*) &temp[26], *(uint64_t*) &temp[36]);

	for (int i = 0; i < inputSize; i += lineLength) {
		uint16_t* count = (uint16_t*) &temp[i];
		printf("===============Count:%"PRIu16"  %i\n", *count, i);
		int iterations = (*count) / 64;
		if ((*count) % 64 > 0) {
			iterations++;
		}
		for (int j = 2; j < iterations * 8; j += 8) {
			printf("%d : %"PRIu64"\n", j, *((uint64_t*) (&temp[i + j])));
		}
		printf("index:%i %"PRIu64", %"PRIu64", %"PRIu64"\n", i,
				*(uint64_t*) &tempFilter[i], *(uint64_t*) &tempFilter[i + 8],
				*(uint64_t*) &tempFilter[i + 16]);
	}

	//		for (int i = 0; i < inputSize; i += lineLength) {
	//			printf("index:%i %"PRIu64", %"PRIu64", %"PRIu64"\n", i,
	//					*(uint64_t*) &tempFilter[i],
	//					*(uint64_t*) &tempFilter[i + 8],
	//					*(uint64_t*) &tempFilter[i + 16]);
	//		}

	//		for (int i = 0; i < inputSize; i++) {
	//			printf("%c", temp[i]);
	//		}

	//		for (int i = 0; i < inputSize; i++) {
	//			printf("%i", (int8_t)tempFilter[i]);
	//		}
	//		printf("\n");

}
