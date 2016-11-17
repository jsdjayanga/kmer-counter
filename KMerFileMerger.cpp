/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   KMerFileMerger.cpp
 * Author: jayanga
 * 
 * Created on November 16, 2016, 2:24 PM
 */

#include <vector>
#include <iostream>
#include <string.h>
#include "KMerFileMerger.h"

KMerFileMerger::KMerFileMerger(list<string> inputFiles, string outputFile, uint64_t kmerLength) {
    
    _kmerStoreLength = kmerLength / 32;
    if (kmerLength % 32 > 0) {
        _kmerStoreLength++;
    }
    _kmerStoreLength *= 8;
    _kmerStoreLength += 4;
    
    for (std::list<string>::iterator it = inputFiles.begin(); it != inputFiles.end(); ++it) {
        _sortedFileList.push_back(new SortedKMerFile(*it, kmerLength));
    }
    
    this->_outputFilename = outputFile;
    
    this->_writeLength = 0;
    this->_writeBufferSize = _kmerStoreLength * 100000;
    this->_writeBuffer = new char[this->_writeBufferSize];
}

KMerFileMerger::KMerFileMerger(const KMerFileMerger& orig) {
}

KMerFileMerger::~KMerFileMerger() {
    for (list<SortedKMerFile*>::iterator it = _processedSortedFileList.begin(); it != _processedSortedFileList.end(); ++it) {
        delete *it;
    }
    delete[] _writeBuffer;
}

void KMerFileMerger::Merge() {
    vector<SortedKMerFile*> lowers;
    while (!_sortedFileList.empty()) {
        std::list<SortedKMerFile*>::iterator it = _sortedFileList.begin();
        SortedKMerFile* sortedKmerFile = *it;

        for (++it; it != _sortedFileList.end(); ++it) {
            SortedKMerFile* currentSortedKmerFile = *it;

            if (CheckLessThan(currentSortedKmerFile->ReadKmer(), sortedKmerFile->ReadKmer(), _kmerStoreLength)) {
                sortedKmerFile = currentSortedKmerFile;
                lowers.clear();
            } else if (CheckEquals(currentSortedKmerFile->ReadKmer(), sortedKmerFile->ReadKmer(), _kmerStoreLength)) {
                lowers.push_back(sortedKmerFile);
                sortedKmerFile = currentSortedKmerFile;
            }
        }
        lowers.push_back(sortedKmerFile);

        SortedKMerFile* firstReader = lowers[0];
        char* firstKmer = firstReader->ReadKmer();
        firstReader->PopKmer();
        for (uint32_t index = 1; index < lowers.size(); index++) {
            SortedKMerFile* kmerReader = lowers[index];
            char* secondKmer = kmerReader->ReadKmer();
            kmerReader->PopKmer();
            *(uint32_t*) (firstKmer + _kmerStoreLength - sizeof (uint32_t)) = *(uint32_t*) (firstKmer + _kmerStoreLength
                    - sizeof (uint32_t)) + *(uint32_t*) (secondKmer + _kmerStoreLength - sizeof (uint32_t));

            if (kmerReader->ReadKmer() == NULL) {
                _processedSortedFileList.push_back(kmerReader);
                _sortedFileList.remove(kmerReader);
            }
        }
        lowers.clear();

        WriteToFile(firstKmer);
        
//        PrintKmer(*(uint64_t*) firstKmer);
//        cout << " " << *(uint64_t*) firstKmer << " " << *(uint32_t*) (firstKmer + 8) << endl;

        if (firstReader->ReadKmer() == NULL) {
            _processedSortedFileList.push_back(firstReader);
            _sortedFileList.remove(firstReader);
        }
    }
    WriteToFile();
}

bool KMerFileMerger::CheckLessThan(char* lhs, char* rhs, uint64_t kmerStoreLength) {
    for (int32_t index = 0; index < kmerStoreLength - sizeof (uint32_t); index += sizeof (uint64_t)) {
        if (*((uint64_t*) (lhs + index)) < *((uint64_t*) (rhs + index))) {
            return true;
        } else if (*((uint64_t*) (lhs + index)) > *((uint64_t*) (rhs + index))) {
            return false;
        }
    }
    return false;
}

bool KMerFileMerger::CheckEquals(char* lhs, char* rhs, uint64_t kmerStoreLength) {
    for (int32_t index = 0; index < kmerStoreLength - sizeof (uint32_t); index += sizeof (uint64_t)) {
        if (*((uint64_t*) (lhs + index)) == *((uint64_t*) (rhs + index))) {
            continue;
        } else {
            return false;
        }
    }
    return true;
}

void KMerFileMerger::WriteToFile(char* record) {
    if (_writeLength + _kmerStoreLength > _writeBufferSize) {
        WriteToFile();
    }
    memcpy(_writeBuffer + _writeLength, record, _kmerStoreLength);
    _writeLength += _kmerStoreLength;
}

void KMerFileMerger::WriteToFile() {
    ofstream out(_outputFilename.c_str(), ios::out | ios::app | ios::binary);
    if (out.is_open()) {
        out.write(_writeBuffer, _writeLength);
        _writeLength = 0;
        out.close();
    }
}
