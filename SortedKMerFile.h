/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   SortedKMerFile.h
 * Author: jayanga
 *
 * Created on November 16, 2016, 12:03 AM
 */

#pragma once

#include <fstream>
#include <stdint.h>
#include <string>

using namespace std;

class SortedKMerFile {
public:
    SortedKMerFile(string filename, uint64_t kmerLenght);
    SortedKMerFile(const SortedKMerFile& orig);
    virtual ~SortedKMerFile();
    char* ReadKmer();
    void PopKmer();
private:
    string _filename;
    uint64_t _kmerLength;
    
    ifstream _fileStream;
    uint64_t _kmerStoreLenght;
    char* _cache;
    uint64_t _cacheSize;
    int64_t _cachePos;
    uint64_t _dataLength;
    bool _accumilated;
    
    void ReadData();
    bool CheckEquals(char* lhs, char* rhs);
    char* GetCurrentKmer();
    char* GetNextKmer();
};


