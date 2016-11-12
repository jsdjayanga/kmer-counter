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

    while (kmerFileReaders.size() > 0) {
        std::list<KMerFileReader*>::iterator it = kmerFileReaders.begin();
        KMerFileReader* reader = *it;

        for (++it; it != kmerFileReaders.end(); ++it) {
            KMerFileReader* currentRader = *it;
            if (checkLessThan((*currentRader).peekKmer(), (*reader).peekKmer(), _kmerLength)) {
                reader = currentRader;
                lowers.clear();
            } else if (checkEquals((*currentRader).peekKmer(), (*reader).peekKmer(), _kmerLength)) {
                lowers.push_back(reader);
                reader = currentRader;
            }
        }

//        printf("===========lowe value=========kmer %"PRIu64", reader=%"PRIu64"\n", ((KMer32*)reader->peekKmer())->kmer[0], reader);
        lowers.push_back(reader);

        KMerFileReader* firstReader = lowers[0];
        char* firstKmer = readWithLocalCount(*firstReader, _kmerLength);
        for (uint32_t index = 1; index < lowers.size(); index++) {
            KMerFileReader* kmerReader = lowers[index];
            char* secondKmer = readWithLocalCount(*kmerReader, _kmerLength);
            ((KMer32*) firstKmer)->count = ((KMer32*) firstKmer)->count + ((KMer32*) secondKmer)->count;
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

bool KMerFileMerger::checkLessThan(char* lhs, char* rhs, uint64_t kmerLength) {
    if (kmerLength <= 32) {
//    	printf("===========checkLessThan=========lhs %"PRIu64", rhs=%"PRIu64"\n", ((KMer32*) lhs)->kmer[0], ((KMer32*) rhs)->kmer[0]);
        return ((KMer32*) lhs)->kmer[0] < ((KMer32*) rhs)->kmer[0];
    }
}

bool KMerFileMerger::checkEquals(char* lhs, char* rhs, uint64_t kmerLength) {
    if (kmerLength <= 32) {
//    	printf("===========checkEquals=========lhs %"PRIu64", rhs=%"PRIu64"\n", ((KMer32*) lhs)->kmer[0], ((KMer32*) rhs)->kmer[0]);
        return ((KMer32*) lhs)->kmer[0] == ((KMer32*) rhs)->kmer[0];
    }
}

char* KMerFileMerger::readWithLocalCount(KMerFileReader& kMerFileReader, uint64_t kmerLength) {
    while (kMerFileReader.peekKmer() != NULL) {
        if (kMerFileReader.peekNextKmer() != NULL && checkEquals(kMerFileReader.peekKmer(), kMerFileReader.peekNextKmer(), kmerLength)) {
            add32Mers(kMerFileReader.peekKmer(), kMerFileReader.peekNextKmer());
            kMerFileReader.popKmer();
        } else {
            return kMerFileReader.popKmer();
        }
    }
    return NULL;
}

void KMerFileMerger::add32Mers(char* first, char* second) {
    ((KMer32*) second)->count = ((KMer32*) second)->count + ((KMer32*) first)->count;
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
