/*
 * KMerFileMerger.cpp
 *
 *  Created on: Nov 10, 2016
 *      Author: jayanga
 */

#include <vector>
#include <cstring>
#include <inttypes.h>
#include "KMerFileMerger.h"
#include "KMerFileHandler.h"
#include "FileDump.h"

KMerFileMerger::KMerFileMerger(string inputDirectory, string outputFilename, uint64_t recordLength,
		uint64_t recordCount, uint64_t kmerLength) {
	this->_inputDirectory = inputDirectory;
	this->_outputFilename = outputFilename;
	this->_recordLength = recordLength;
	this->_recordCount = recordCount;
	this->_kmerLength = kmerLength;
	this->_writeBufferSize = recordLength * recordCount * 100;
	this->_writeBuffer = new char[this->_writeBufferSize];
	this->_writeLength = 0;
	this->_kMerFileHandler = new KMerFileHandler(this->_inputDirectory, this->_recordLength, this->_recordCount);

}

KMerFileMerger::~KMerFileMerger() {
	// TODO Auto-generated destructor stub
}

void KMerFileMerger::merge() {

	list<KMerFileReader*>& kmerFileReaders = this->_kMerFileHandler->getKmerFileList();

	vector<KMerFileReader*> lowers;
	uint32_t fileCount = kmerFileReaders.size();

	uint64_t kmerStoreSize = _kmerLength / 32;
	if (_kmerLength % 32 > 0) {
		kmerStoreSize++;
	}
	kmerStoreSize *= 8;
	kmerStoreSize += 4;

	while (kmerFileReaders.size() > 0) {
		std::list<KMerFileReader*>::iterator it = kmerFileReaders.begin();
		KMerFileReader* reader = *it;

		for (++it; it != kmerFileReaders.end(); ++it) {
			KMerFileReader* currentRader = *it;

//			printf("================= currentReaderKM=%"PRIu64", ReaderKM=%"PRIu64", lessThan=%d\n",
//					*(uint64_t*) currentRader->peekKmer(), *(uint64_t*) reader->peekKmer(),
//					*(uint64_t*) currentRader->peekKmer() < *(uint64_t*) reader->peekKmer());

			if (checkLessThan((*currentRader).peekKmer(), (*reader).peekKmer(), kmerStoreSize)) {
//				printf("===================Less than\n");
				reader = currentRader;
				lowers.clear();
			} else if (checkEquals((*currentRader).peekKmer(), (*reader).peekKmer(), kmerStoreSize)) {
//				printf("===================Equals========\n");
				lowers.push_back(reader);
				reader = currentRader;
			}
		}

//        printf("===========lowe value=========kmer %"PRIu64", reader=%"PRIu64"\n", ((KMer32*)reader->peekKmer())->kmer[0], reader);
		lowers.push_back(reader);

//		printf("====================Selection Done=======Lowers Count = %i\n", lowers.size());

		KMerFileReader* firstReader = lowers[0];
		char* firstKmer = readWithLocalCount(*firstReader, kmerStoreSize);
		for (uint32_t index = 1; index < lowers.size(); index++) {
			KMerFileReader* kmerReader = lowers[index];
			char* secondKmer = readWithLocalCount(*kmerReader, kmerStoreSize);
			//((KMer32*) firstKmer)->count = ((KMer32*) firstKmer)->count + ((KMer32*) secondKmer)->count;

			*(uint32_t*) (firstKmer + kmerStoreSize - sizeof(uint32_t)) = *(uint32_t*) (firstKmer + kmerStoreSize
					- sizeof(uint32_t)) + *(uint32_t*) (secondKmer + kmerStoreSize - sizeof(uint32_t));

			kmerReader->popKmer();
			if (kmerReader->peekKmer() == NULL) {
				kmerFileReaders.remove(kmerReader);
			}
		}
		lowers.clear();

		writeToFile(firstKmer);

		if (firstReader->peekKmer() == NULL) {
			kmerFileReaders.remove(firstReader);
		}

		// TODO : Remove Test Codes ==============================
//        FileDump* fileDump = new FileDump("/tmp/mydata.dump");
//        fileDump->dumpToConsole2(_writeBuffer, 32, _writeLength);
//        printf("==========================================================\n");
		//========================================================

	}

	writeToFile();
}

bool KMerFileMerger::checkLessThan(char* lhs, char* rhs, uint64_t kmerStoreSize) {
	for (int32_t index = 0; index < kmerStoreSize - sizeof(uint32_t); index += sizeof(uint64_t)) {
//		printf("================= lhs=%"PRIu64", rhs=%"PRIu64"\n", *(uint64_t*) (lhs + index),
//				*(uint64_t*) (rhs + index));
		if (*((uint64_t*) (lhs + index)) < *((uint64_t*) (rhs + index))) {
			return true;
		} else if (*((uint64_t*) (lhs + index)) > *((uint64_t*) (rhs + index))) {
			return false;
		}
	}
	return false;
}

bool KMerFileMerger::checkEquals(char* lhs, char* rhs, uint64_t kmerStoreSize) {
	for (int32_t index = 0; index < kmerStoreSize - sizeof(uint32_t); index += sizeof(uint64_t)) {
		if (*((uint64_t*) (lhs + index)) == *((uint64_t*) (rhs + index))) {
			continue;
		} else {
			return false;
		}
	}
	return true;
}

char* KMerFileMerger::readWithLocalCount(KMerFileReader& kMerFileReader, uint64_t kmerStoreSize) {
	while (kMerFileReader.peekKmer() != NULL) {
		if (kMerFileReader.peekNextKmer() != NULL
				&& checkEquals(kMerFileReader.peekKmer(), kMerFileReader.peekNextKmer(), kmerStoreSize)) {
			//add32Mers(kMerFileReader.peekKmer(), kMerFileReader.peekNextKmer());
			*(uint32_t*) (kMerFileReader.peekNextKmer() + kmerStoreSize - sizeof(uint32_t)) =
					*(uint32_t*) (kMerFileReader.peekNextKmer() + kmerStoreSize - sizeof(uint32_t))
							+ *(uint32_t*) (kMerFileReader.peekKmer() + kmerStoreSize - sizeof(uint32_t));
			kMerFileReader.popKmer();
		} else {
			return kMerFileReader.popKmer();
		}
	}
	return NULL;
}

void KMerFileMerger::writeToFile(char* record) {
	if (_writeLength + _recordLength > _writeBufferSize) {
		writeToFile();
	}
	memcpy(_writeBuffer + _writeLength, record, _recordLength);
	_writeLength += _recordLength;
}

void KMerFileMerger::writeToFile() {

	// TODO : Remove Test Codes ==============================
//    FileDump* fileDump = new FileDump("/tmp/mydata.dump");
//    fileDump->dumpToConsole2(_writeBuffer, 32, _writeLength);
	//========================================================

	ofstream out(_outputFilename.c_str(), ios::out | ios::app | ios::binary);
	if (out.is_open()) {
		out.write(_writeBuffer, _writeLength);
		_writeLength = 0;
		out.close();
	}
}
