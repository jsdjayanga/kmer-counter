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

#include <stdint.h>

using namespace std;



struct KMer32
{
    uint64_t kmer[1];
    uint32_t count;
};

struct KMer64
{
    uint64_t kmer[2];
    uint32_t count;
};

struct KMer96
{
    uint64_t kmer[3];
    uint32_t count;
};

struct KMer128
{
    uint64_t kmer[4];
    uint32_t count;
};

int64_t processKMers(const char* input, int64_t kmerLength, int64_t inputSize,
		int64_t lineLength);

void printBitEncodedResult(char* d_input, char* d_filter, uint64_t inputSize,uint64_t lineLength);
void printKmerResult(char* d_output, uint64_t outputSize, uint64_t kmerLength);
void dumpKmersWithLengthToConsole(char* data, int64_t lineLength, int64_t outputLength, uint64_t kmerLenght);
void printDNABase(uint64_t value);

