/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Options.cpp
 * Author: jayanga
 * 
 * Created on October 19, 2016, 8:05 PM
 */

#include "Options.h"

Options::Options() {
    this->_numberOfThreads = 0;
    this->_chunkSize = 10000;
    this->_gpuMemoryLimit = 10000000;
    this->_kmerLength = 32;
}

Options::Options(const Options& orig) {
}

Options::~Options() {
}

void Options::SetInputFileDirectory(string inputFileDirectory) {
    this->_inputFileDirectory = inputFileDirectory;
}

string Options::GetInputFileDirectory() const {
    return _inputFileDirectory;
}

void Options::SetChunkSize(int64_t _chunkSize) {
    this->_chunkSize = _chunkSize;
}

int64_t Options::GetChunkSize() const {
    return _chunkSize;
}

void Options::SetNumberOfThreads(int64_t _numberOfThreads) {
    this->_numberOfThreads = _numberOfThreads;
}

int64_t Options::GetNumberOfThreads() const {
    return _numberOfThreads;
}

void Options::SetGpuMemoryLimit(int64_t _gpuMemoryLimit) {
    this->_gpuMemoryLimit = _gpuMemoryLimit;
}

int64_t Options::GetGpuMemoryLimit() const {
    return _gpuMemoryLimit;
}

void Options::SetKmerLength(int64_t _kmerLength) {
    this->_kmerLength = _kmerLength;
}

int64_t Options::GetKmerLength() const {
    return _kmerLength;
}


