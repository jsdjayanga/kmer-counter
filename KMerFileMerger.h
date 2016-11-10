/*
 * KMerFileMerger.h
 *
 *  Created on: Nov 10, 2016
 *      Author: jayanga
 */

#pragma once

#include <list>
#include "KMerFileReader.h"
#include "KMerFileHandler.h"
#include "KMerSizes.h"

using namespace std;

class KMerFileMerger {
public:
	KMerFileMerger(string inputDirectory, string outputFilename, uint64_t recordLength, uint64_t recordCount,
			uint64_t kmerLength);
	virtual ~KMerFileMerger();
	void merge();
private:
	string _inputDirectory;
	string _outputFilename;
	uint64_t _recordLength;
	uint64_t _recordCount;
	uint64_t _kmerLength;

	KMerFileHandler* _kMerFileHandler;
	char* _writeBuffer;
	uint64_t _writeBufferSize;
	uint64_t _writeLength;

	bool checkLessThan(char* lhs, char* rhs, uint64_t kmerLength);
	bool checkEquals(char* lhs, char* rhs, uint64_t kmerLength);
	char* readWithLocalCount(KMerFileReader& kMerFileReader, uint64_t kmerLength);
	void writeToFile(char* record);
	void writeToFile();

	void add32Mers(char* first, char* second);
};

