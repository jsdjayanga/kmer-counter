/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   GPUHandler.h
 * Author: jayanga
 *
 * Created on October 21, 2016, 10:22 PM
 */

#pragma once
#pragma pack(1)

#include <stdlib.h>
#include <stdint.h>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "FileDump.h"
#include "KMerSizes.h"

using namespace std;

#define cudaErrorCheck(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
	if (code != cudaSuccess) {
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort)
			exit(code);
	}
}

int64_t processKMers(const char* input, int64_t kmerLength, int64_t inputSize, int64_t lineLength, uint32_t readId,
		FileDump& fileDump);

void printBitEncodedResult(char* d_input, char* d_filter, uint64_t inputSize, uint64_t lineLength);
void printKmerResult(char* d_output, uint64_t outputSize, uint64_t kmerLength);
void dumpKmersWithLengthToConsole(char* data, int64_t lineLength, int64_t outputLength, uint64_t kmerLenght);
void dumpKmersWithLengthToConsoleHost(char* data, int64_t lineLength, int64_t outputLength, uint64_t kmerLenght);
void printDNABase(uint64_t value);

