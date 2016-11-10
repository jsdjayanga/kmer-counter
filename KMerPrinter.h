/*
 * KMerPrinter.h
 *
 *  Created on: Nov 10, 2016
 *      Author: jayanga
 */

#pragma once

#include <string>
#include <stdint.h>

using namespace std;

class KMerPrinter {
public:
	KMerPrinter(string filename, uint64_t kmerlength);
	virtual ~KMerPrinter();
	void print();
private:
	string _filename;
	uint64_t _kmerlength;
	uint64_t _kmerStorelength;
	uint64_t _chunkSize;
	char* _data;
	uint64_t _availableData;

	void printKmer(uint64_t value);
};

