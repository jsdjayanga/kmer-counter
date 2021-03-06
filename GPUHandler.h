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
#include <list>
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

struct GPUStream {
	GPUStream(uint32_t id, char* h_output, char* d_input, char* d_output, char* d_filter) {
		_id = id;
		_h_output = h_output;
		_d_input = d_input;
		_d_output = d_output;
		_d_filter = d_filter;

		_kmer_db.clear();
		_kmer_db_line_length = 250 * 1024 * 1024;
		_kmer_db.push_front(new char[_kmer_db_line_length]);
		_kmer_db_line_index = 0;
	}
	cudaStream_t stream;
	char* _h_output;
	char* _d_input;
	char* _d_output;
	char* _d_filter;
	uint32_t _id;

	list<char*> _kmer_db;
	uint64_t _kmer_db_line_length;
	uint64_t _kmer_db_line_index;
};

GPUStream** PrepareGPU(uint32_t streamCount, uint64_t inputSize, uint64_t lineLength, int64_t kmerLength);
void FreeGPU(GPUStream** streams, uint32_t streamCount);

int64_t processKMers(GPUStream* gpuStream, const char* input, int64_t kmerLength, int64_t inputSize, int64_t lineLength, uint32_t readId,
		FileDump& fileDump);

void printBitEncodedResult(char* d_input, char* d_filter, uint64_t inputSize, uint64_t lineLength);
void printKmerResult(char* d_output, uint64_t outputSize, uint64_t kmerLength);
void dumpKmersWithLengthToConsole(char* data, int64_t lineLength, int64_t outputLength, uint64_t kmerLenght);
void dumpKmersWithLengthToConsoleHost(char* data, int64_t lineLength, int64_t outputLength, uint64_t kmerLenght);
void printDNABase(uint64_t value);

