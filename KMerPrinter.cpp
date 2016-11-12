/*
 * KMerPrinter.cpp
 *
 *  Created on: Nov 10, 2016
 *      Author: jayanga
 */

#include <string.h>
#include <fstream>
#include <inttypes.h>
#include "KMerPrinter.h"

KMerPrinter::KMerPrinter(string inputFilename, string outputFilename, uint64_t kmerlength) {
	this->_inputFilename = inputFilename;
	this->_outputFilename = outputFilename;
	this->_kmerlength = kmerlength;
	this->_availableData = 0;

	int entryLength = kmerlength / 32;
	if (kmerlength % 32 > 0) {
		entryLength++;
	}
	entryLength *= 8;
	entryLength += 4;
	this->_kmerStorelength = entryLength;
	this->_chunkSize = this->_kmerStorelength * 10000;
	//this->_chunkSize = 50000000;
	this->_data = new char[this->_chunkSize];
}

KMerPrinter::~KMerPrinter() {
	// TODO Auto-generated destructor stub
}

void KMerPrinter::print() {
	ifstream fileStream;
	fileStream.open(_inputFilename.c_str(), ios::binary);
	if (fileStream.is_open()) {
		while (!fileStream.eof()) {
			memset(_data, 0, _chunkSize);
			fileStream.read(_data, _chunkSize);
			_availableData = fileStream.gcount();
//			printf("=================================================reading a chunk form file ==%"PRIu64"\n", _availableData);
//
//			for (uint64_t ind = 0; ind < _availableData; ind += 12) {
//				printf("====print==%"PRIu64", %u\n", *(uint64_t*) (_data + ind), *(uint32_t*) (_data + ind + 8));
//			}

			for (int64_t index = 0; index < _availableData; index += _kmerStorelength) {

				int kmerBytes = _kmerStorelength - 4;
				int i = index;
				for (; i < index + kmerBytes; i += 8) {
					uint64_t value = 0;
					memcpy(&value, &_data[i], sizeof(uint64_t));
					printKmer(value);
					printf(" %"PRIu64" ", value);
				}

				uint32_t count = 0;
				memcpy(&count, &_data[i], sizeof(uint32_t));
				printf(" %u\n", count);
			}
		}
	}
}

void KMerPrinter::printKmer(uint64_t value) {
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

