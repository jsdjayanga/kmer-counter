/*
 * KMerFileReader.h
 *
 *  Created on: Nov 10, 2016
 *      Author: jayanga
 */

#pragma once

#include <string>
#include <stdint.h>
#include <fstream>

using namespace std;

class KMerFileReader {
public:
	KMerFileReader(string filename, uint64_t fileSize, uint64_t recordLength, uint64_t recordCount);
	virtual ~KMerFileReader();
	char* peekKmer();
	char* peekNextKmer();
	char* popKmer();
private:
	ifstream _fileStream;
	string _filename;
	uint64_t _fileSize;
	uint64_t _recordLength;
	uint64_t _recordCount;

	char* _cache;
	uint64_t _cachePos;
	uint64_t _cacheSize;
	uint64_t _availableLength;

	void readData(uint64_t cacheStartPos, uint64_t size);
};

