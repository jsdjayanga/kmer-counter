#include <stdio.h>
#include <inttypes.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include "GPUHandler.h"

__global__ void bitEncode(char* input, char* filter, int64_t lineLength, int64_t upperBound) {
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
			int64_t readValueLocation = ((((i - index) / 32) - 1) * sizeof(int64_t)) + sizeof(int16_t) + index;
			//int64_t filterLocation = ((((i - index) / 32) - 1) * sizeof (int64_t));
			memcpy(&input[readValueLocation], &readValue, sizeof(int64_t));
			//memcpy(&filter[filterLocation], &filterValue, sizeof (int64_t));
			readValue = 0;

			//filterValue = 0;
		}

		if (i - index > 0 && (i - index) % 64 == 0) {
			int64_t filterLocation = ((((i - index) / 64) - 1) * sizeof(int64_t)) + index;
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

	//printf("readValueLength==============:%"PRIu16"\n", readValueLength);
	memcpy(&input[index], &readValueLength, sizeof(uint16_t));

	if (i > 0 && (i - index) % 64 > 0) {
		uint8_t shiftingReadValue = (32 - ((i - index) % 32)) * 2;
		readValue <<= shiftingReadValue;

		int64_t readValueLocation = ((((i - index) / 32)) * sizeof(int64_t)) + sizeof(int16_t) + index;
		memcpy(&input[readValueLocation], &readValue, sizeof(int64_t));
		readValue = 0;
		readValueLength = 0;

		uint8_t shiftingFilterValue = 64 - ((i - index) % 64);
		filterValue <<= shiftingFilterValue;
		int64_t filterLocation = ((((i - index) / 64)) * sizeof(int64_t)) + index;
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

__global__ void extractKMers(char* input, char* bitFilter, char*output, uint64_t sectionLength, int64_t kmerLength,
		int64_t upperBound, int64_t lineLength) {
	uint64_t index = (blockIdx.x * blockDim.x + threadIdx.x) * sectionLength;
	uint64_t filterIndex = (blockIdx.x * blockDim.x + threadIdx.x) * lineLength;

	if (filterIndex > upperBound - lineLength) {
		return;
	}

//printf("================================extract=========================================:index=%"PRIu64"\n", filterIndex);

//	uint64_t* f1 = (uint64_t*) &temp[0];
//	uint64_t* f2 = (uint64_t*) &temp[8];
//	uint64_t* f3 = (uint64_t*) &temp[16];
//	cout << *f1 << "|" << *f2 << "|" << *f3 << endl;

	uint16_t i1 = *(uint16_t*) &input[filterIndex];
	uint16_t filterLength = i1 / 2;

	char* encodedInput = &input[filterIndex + sizeof(uint16_t)];

//	cout << "=========:filter length=" << filterLengrh << endl;
	int64_t filterReadLength = 0;
//printf("=======================b4 for loop, filterIndex=%"PRIu64", filterLength=%"PRIu64"\n", filterIndex, filterLength);
	uint64_t outputIndex = 0;
	//uint64_t i = filterIndex;
	for (uint64_t i = 0; i < filterLength; i++) {

//printf("=======================inside for loopindex=%"PRIu64"\n", filterIndex);
		uint64_t filter = 0;
		if (i == 0) {
			memcpy(&filter, &bitFilter[filterIndex], sizeof(uint64_t));
			//filter = (uint64_t*) &bitFilter[filterIndex];
		} else {
			memcpy(&filter, &bitFilter[((i / 64) * 8) + filterIndex], sizeof(uint64_t));
			//filter = (uint64_t*) &bitFilter[(((filterIndex + i) / 64) * 8)];
		}
//		cout << "=====filter::::" << *filter << "|" << i << "|"
//				<< checkBit(*filter, (uint8_t) i % 64) << endl;

		if (!checkBit(filter, (uint8_t) i % 64)) {
			filterReadLength++;

			if (filterReadLength >= kmerLength) {
//				cout << "valid kmer" << endl;
//printf("================================validkmeri=%"PRIu64", index=%"PRIu64"\n", i, filterIndex);
				uint64_t firstByte = i - kmerLength + 1; // +1 is needed as 'i' starts with index 0
				int64_t shifting = ((firstByte % 32)) * 2;
				int64_t rightShifting = (32 - (kmerLength % 32)) * 2;

				int64_t firstByteToReadEncodedInput = ((firstByte / 32) * 8);
//printf("==========================================================firstByteToReadEncodedInput %"PRIu64", index=%"PRIu64"\n", firstByteToReadEncodedInput, filterIndex);
				uint16_t kmerByteLength = kmerLength / 4;
				uint64_t kmerByteStoreLength = kmerLength / 4;
				if (kmerLength % 4 > 0) {
					kmerByteLength++;
					kmerByteStoreLength += sizeof(uint64_t);
				}

				uint64_t readValue1 = 0;
				uint64_t readValue2 = 0;

				for (int64_t x = firstByteToReadEncodedInput; x < firstByteToReadEncodedInput + kmerByteLength; x +=
						sizeof(uint64_t)) {

					//memcpy(&readValue, (char*) &encodedInput[x], sizeof (uint64_t));
					readValue1 = read64bits(encodedInput, x);
//printf("==========================================================actualReadingIndex %"PRIu64", index=%"PRIu64", readValue1Ori=%"PRIu64", shifting=%"PRIu64"\n", x, filterIndex, readValue1, shifting);
					if (shifting > 0) {
						readValue1 <<= shifting;

						if ((x + sizeof(uint64_t)) * 4 < filterLength) { // To avoid reading beyond the limit of bit encoded input
							readValue2 = read64bits(encodedInput, x + sizeof(uint64_t));
							readValue2 >>= (64 - shifting);
						}
						readValue1 |= readValue2;
					}

//					cout << "=============================AAAAAA:" << x << "|"
//							<< readValue1 << "|" << index + outputIndex << endl;

					if (x + sizeof(uint64_t) > firstByteToReadEncodedInput + kmerByteLength) {
						readValue1 >>= rightShifting;
						readValue1 <<= rightShifting;
					}

					memcpy(&output[index + outputIndex], &readValue1, sizeof(uint64_t));
//printf("==========================================================readValue1 %"PRIu64", actualOutputIndex=%"PRIu64"\n", readValue1, index + outputIndex);
					outputIndex += sizeof(uint64_t);

					readValue1 = 0;
					readValue2 = 0;
				}
				uint32_t count = 1;
				memcpy(&output[index + outputIndex], &count, sizeof(uint32_t));
				outputIndex += sizeof(uint32_t);

				filterReadLength--;
			}

		} else {
			filterReadLength = 0;
		}
	}
}

uint64_t calculateOutputSize(int64_t inputSize, int64_t lineLength, int64_t kmerLength) {
	uint64_t records = inputSize / lineLength;
	uint64_t kmerCount = lineLength - kmerLength + 1;
	uint64_t kmerStoreSize = kmerLength / 32;
	if (kmerLength % 32 > 0) {
		kmerStoreSize++;
	}
	kmerStoreSize *= 8;
	kmerStoreSize += 4;
	return kmerCount * kmerStoreSize * records;
}

class KMer32Comparator {
public:
	__device__ __host__
	bool operator()(const KMer32 &c1, const KMer32 &c2) {
		return c1.kmer[0] < c2.kmer[0];
	}
};

class KMer64Comparator {
public:
	__device__ __host__
	bool operator()(const KMer64 &c1, const KMer64 &c2) {
		for (int32_t index = 0; index < (sizeof(c1.kmer) / sizeof(*c1.kmer)); index++) {
			if (c1.kmer[index] < c2.kmer[index]) {
				return true;
			} else if (c1.kmer[index] > c2.kmer[index]) {
				return false;
			}
		}
		return false;
	}
};

class KMer96Comparator {
public:
	__device__ __host__
	bool operator()(const KMer96 &c1, const KMer96 &c2) {
		for (int32_t index = 0; index < (sizeof(c1.kmer) / sizeof(*c1.kmer)); index++) {
			if (c1.kmer[index] < c2.kmer[index]) {
				return true;
			} else if (c1.kmer[index] > c2.kmer[index]) {
				return false;
			}
		}
		return false;
	}
};

class KMer128Comparator {
public:
	__device__ __host__
	bool operator()(const KMer128 &c1, const KMer128 &c2) {
		for (int32_t index = 0; index < (sizeof(c1.kmer) / sizeof(*c1.kmer)); index++) {
			if (c1.kmer[index] < c2.kmer[index]) {
				return true;
			} else if (c1.kmer[index] > c2.kmer[index]) {
				return false;
			}
		}
		return false;
	}
};

void sortKmers(char* d_output, uint64_t kmerLength, uint64_t outputSize) {
	uint64_t kmerStoreSize = kmerLength / 32;
	if (kmerLength % 32 > 0) {
		kmerStoreSize++;
	}
	kmerStoreSize *= 8;
	kmerStoreSize += 4;

	if (kmerLength <= 32) {
		printf("invoking ==== KMer32Comparator");
		thrust::device_ptr<KMer32> d_Kmers((KMer32*) d_output);
		thrust::sort(d_Kmers, d_Kmers + (outputSize / kmerStoreSize), KMer32Comparator());
	} else if (kmerLength <= 64) {
		printf("invoking ==== KMer64Comparator");
		thrust::device_ptr<KMer64> d_Kmers((KMer64*) d_output);
		thrust::sort(d_Kmers, d_Kmers + (outputSize / kmerStoreSize), KMer64Comparator());
	} else if (kmerLength <= 96) {
		printf("invoking ==== KMer96Comparator");
		thrust::device_ptr<KMer96> d_Kmers((KMer96*) d_output);
		thrust::sort(d_Kmers, d_Kmers + (outputSize / kmerStoreSize), KMer96Comparator());
	} else if (kmerLength <= 128) {
		printf("invoking ==== KMer128Comparator");
		thrust::device_ptr<KMer128> d_Kmers((KMer128*) d_output);
		thrust::sort(d_Kmers, d_Kmers + (outputSize / kmerStoreSize), KMer128Comparator());
	} else {
		printf("Sorting is not supported for kmers with length higher than %"PRIu64"\n", kmerLength);
	}
}

int64_t processKMers(const char* input, int64_t kmerLength, int64_t inputSize, int64_t lineLength, uint32_t readId,
		FileDump& fileDump) {
	printf("Processing k-mers klen=%"PRIu64", inSize=%"PRIu64","
	" liLen=%"PRIu64"\n", kmerLength, inputSize, lineLength);

	char* h_output;

	char* d_input;
	char* d_output;
	char* d_filter;

	uint64_t outputSize = calculateOutputSize(inputSize, lineLength, kmerLength);
	printf("Output size =============%"PRIu64"\n", outputSize);

	h_output = (char *) malloc(outputSize);
	memset(h_output, 0, outputSize);

	cudaErrorCheck(cudaMalloc((void ** ) &d_input, inputSize));
	cudaErrorCheck(cudaMalloc((void ** ) &d_output, outputSize));
	cudaErrorCheck(cudaMalloc((void ** ) &d_filter, inputSize));

	cudaErrorCheck(cudaMemcpy(d_input, input, inputSize, cudaMemcpyHostToDevice));
	cudaErrorCheck(cudaMemset(d_output, 0, outputSize));
	cudaErrorCheck(cudaMemset(d_filter, 0, inputSize));

	int32_t threadCount = 512;
	int32_t blockCount = 512;
	int32_t totalThread = blockCount * threadCount;
	int32_t count = inputSize / lineLength / totalThread;
	if ((inputSize / lineLength) % threadCount > 0) {
		count++;
	}

	for (int32_t ite = 0; ite < count; ite++) {
		bitEncode<<<blockCount, threadCount>>>(&d_input[threadCount * lineLength * ite],
				&d_filter[threadCount * lineLength * ite],
				lineLength,
				inputSize);
		cudaErrorCheck(cudaPeekAtLastError());
		cudaErrorCheck(cudaDeviceSynchronize());

		extractKMers<<<blockCount, threadCount>>>(
				&d_input[threadCount * lineLength * ite],
				&d_filter[threadCount * lineLength * ite],
				&d_output[threadCount * outputSize / (inputSize / lineLength) * ite],
				outputSize / (inputSize / lineLength),
				kmerLength,
				inputSize,
				lineLength);
		cudaErrorCheck(cudaPeekAtLastError());
		cudaErrorCheck(cudaDeviceSynchronize());
	}

	printf("Before Sort lineLength=%"PRIu64", outputSize=%"PRIu64", kmerLength=%"PRIu64"\n", lineLength, outputSize,
			kmerLength);
	//dumpKmersWithLengthToConsole(d_output, lineLength, outputSize, kmerLength);

	// Sort step
	sortKmers(d_output, kmerLength, outputSize);
	cudaErrorCheck(cudaPeekAtLastError());
	cudaErrorCheck(cudaDeviceSynchronize());

	cudaErrorCheck(cudaMemcpy(h_output, d_output, outputSize, cudaMemcpyDeviceToHost));
	fileDump.dumpKmersToFile(readId, h_output, outputSize);
	//printBitEncodedResult(d_input, d_filter, inputSize, lineLength);

	//printKmerResult(d_output, outputSize, kmerLength);
	printf("After Sort\n");
	dumpKmersWithLengthToConsoleHost(h_output, lineLength, outputSize, kmerLength);

	cudaFree(d_input);
	cudaFree(d_output);
	cudaFree(d_filter);
	cudaDeviceReset();

	return 0;
}

void printBitEncodedResult(char* d_input, char* d_filter, uint64_t inputSize, uint64_t lineLength) {
	char* temp = new char[inputSize];
	memset(temp, 0, inputSize);

	char* tempFilter = new char[inputSize];
	memset(tempFilter, 0, inputSize);

	cudaMemcpy(temp, d_input, inputSize, cudaMemcpyDeviceToHost);
	cudaMemcpy(tempFilter, d_filter, inputSize, cudaMemcpyDeviceToHost);

	printf("%"PRIu16" : %"PRIu64", %"PRIu64", %"PRIu64", %"PRIu64", %"PRIu64"\n", *(uint16_t*) &temp[0],
			*(uint64_t*) &temp[2], *(uint64_t*) &temp[10], *(uint64_t*) &temp[18], *(uint64_t*) &temp[26],
			*(uint64_t*) &temp[36]);

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
		printf("index:%i %"PRIu64", %"PRIu64", %"PRIu64"\n", i, *(uint64_t*) &tempFilter[i],
				*(uint64_t*) &tempFilter[i + 8], *(uint64_t*) &tempFilter[i + 16]);
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

void printKmerResult(char* d_output, uint64_t outputSize, uint64_t kmerLength) {
	char* temp = new char[outputSize];
	memset(temp, 0, outputSize);

	cudaMemcpy(temp, d_output, outputSize, cudaMemcpyDeviceToHost);

	uint64_t kmerByteLength = kmerLength / 32;
	if (kmerLength % 32 > 0) {
		kmerByteLength += 1;
	}
	kmerByteLength *= 8;

	for (int i = 0; i < outputSize; i += kmerByteLength + 4) {
		int j = 0;
		for (; j < kmerByteLength; j += 8) {
			printf("%"PRIu64" ", *(uint64_t*) &temp[i + j]);
		}
		printf("%"PRIu32"\n", *(uint32_t*) &temp[i + j]);
	}

}

void dumpKmersWithLengthToConsoleHost(char* data, int64_t lineLength, int64_t outputSize, uint64_t kmerLenght) {
	int entryLength = kmerLenght / 32;
	if (kmerLenght % 32 > 0) {
		entryLength++;
	}
	entryLength *= 8;
	entryLength += 4;

	for (int64_t index = 0; index < outputSize; index += entryLength) {
		int kmerBytes = entryLength - 4;
		int i = index;
		for (; i < index + kmerBytes; i += 8) {
			uint64_t value = 0;
			memcpy(&value, &data[i], sizeof(uint64_t));
			printDNABase(value);
			printf(" %" PRIu64 " ", value);
		}

		uint32_t count = 0;
		memcpy(&count, &data[i], sizeof(uint32_t));
		printf(" %u\n", count);
	}
}

void dumpKmersWithLengthToConsole(char* d_data, int64_t lineLength, int64_t outputSize, uint64_t kmerLenght) {

	char* data = new char[outputSize];
	memset(data, 0, outputSize);

	cudaMemcpy(data, d_data, outputSize, cudaMemcpyDeviceToHost);

	int entryLength = kmerLenght / 32;
	if (kmerLenght % 32 > 0) {
		entryLength++;
	}
	entryLength *= 8;
	entryLength += 4;

	for (int64_t index = 0; index < outputSize; index += entryLength) {
		int kmerBytes = entryLength - 4;
		int i = index;
		for (; i < index + kmerBytes; i += 8) {
			uint64_t value = 0;
			memcpy(&value, &data[i], sizeof(uint64_t));
			printDNABase(value);
			printf(" %" PRIu64 " ", value);
		}

		uint32_t count = 0;
		memcpy(&count, &data[i], sizeof(uint32_t));
		printf(" %u\n", count);
	}
}

void printDNABase(uint64_t value) {
	for (int64_t index = 0; index < 64; index += 2) {
		uint64_t temp = value << index;
		temp >>= 62;

		switch (temp) {
		case 0:
			printf("A");
			continue;
		case 1:
			printf("C");
			continue;
		case 2:
			printf("G");
			continue;
		case 3:
			printf("T");
			continue;
		default:
			printf("-");
			continue;
		}
	}
}
