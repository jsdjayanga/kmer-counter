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
	const string& getOutputFile() const;
	void setOutputFile(const string& outputFile);
	const string& getTempFileLocation() const;
	void setTempFileLocation(const string& tempFileLocation);
	uint32_t getNoOfMergersAtOnce() const;
	void setNoOfMergersAtOnce(uint32_t noOfMergersAtOnce);
	uint32_t getNoOfMergeThreads() const;
	void setNoOfMergeThreads(uint32_t noOfMergeThreads);

private:
	string _inputFileDirectory;
	int64_t _numberOfThreads;
	int64_t _chunkSize;
	int64_t _gpuMemoryLimit;
	int64_t _kmerLength;

	string _tempFileLocation;
	string _outputFile;

	uint32_t _noOfMergersAtOnce;
	uint32_t _noOfMergeThreads;
};

