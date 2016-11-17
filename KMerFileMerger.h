/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   KMerFileMerger.h
 * Author: jayanga
 *
 * Created on November 16, 2016, 2:24 PM
 */

#pragma once

#include <string>
#include <list>
#include <stdint.h>
#include "SortedKMerFile.h"

using namespace std;

class KMerFileMerger {
public:
    KMerFileMerger(list<string> inputFiles, string outputFile, uint64_t kmerLength);
    KMerFileMerger(const KMerFileMerger& orig);
    virtual ~KMerFileMerger();

    void Merge();
private:
    bool CheckLessThan(char* lhs, char* rhs, uint64_t kmerStoreLength);
    bool CheckEquals(char* lhs, char* rhs, uint64_t kmerStoreLength);
    void WriteToFile(char* record);
    void WriteToFile();

    void PrintKmer(uint64_t value);
    
    uint64_t _kmerStoreLength;
    char* _writeBuffer;
    uint64_t _writeBufferSize;
    uint64_t _writeLength;
    
    list<SortedKMerFile*> _sortedFileList;
    list<SortedKMerFile*> _processedSortedFileList;
    string _outputFilename;
};


