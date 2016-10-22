/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Options.h
 * Author: jayanga
 *
 * Created on October 19, 2016, 8:05 PM
 */

#pragma once

#include <string>
#include <stdint.h>

using namespace std;

class Options {
public:
    Options();
    Options(const Options& orig);
    virtual ~Options();
    void SetInputFileDirectory(string inputFileDirectory);
    string GetInputFileDirectory() const;
    void SetChunkSize(int64_t _chunkSize);
    int64_t GetChunkSize() const;
    void SetNumberOfThreads(int64_t _numberOfThreads);
    int64_t GetNumberOfThreads() const;
    void SetGpuMemoryLimit(int64_t _gpuMemoryLimit);
    int64_t GetGpuMemoryLimit() const;
    void SetKmerLength(int64_t _kmerLength);
    int64_t GetKmerLength() const;
    
private:
    string  _inputFileDirectory;
    int64_t _numberOfThreads;
    int64_t _chunkSize;
    int64_t _gpuMemoryLimit;
    int64_t _kmerLength;
};


