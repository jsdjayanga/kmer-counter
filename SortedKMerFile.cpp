/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   SortedKMerFile.cpp
 * Author: jayanga
 * 
 * Created on November 16, 2016, 12:03 AM
 */

#include <string.h>
#include <iostream>
#include "SortedKMerFile.h"

SortedKMerFile::SortedKMerFile(string filename, uint64_t kmerLenght) {
    this->_filename = filename;
    this->_kmerLength = kmerLenght;

    this->_kmerStoreLenght = this->_kmerLength / 32;
    if (this->_kmerLength % 32 > 0) {
        _kmerStoreLenght++;
    }
    this->_kmerStoreLenght *= 8;
    this->_kmerStoreLenght += 4;

    this->_cacheSize = this->_kmerStoreLenght * 10000;
    this->_cache = new char[this->_cacheSize];
    memset(this->_cache, 0, this->_cacheSize);

    this->_fileStream.open(filename.c_str());
    this->_fileStream.seekg(0, ios_base::beg);
    this->_cachePos = 0;
    this->_dataLength = 0;

    this->_accumilated = false;

    ReadData();
}

SortedKMerFile::SortedKMerFile(const SortedKMerFile& orig) {
}

SortedKMerFile::~SortedKMerFile() {
    delete[] _cache;
}

void SortedKMerFile::ReadData() {
    memmove(_cache, _cache + _cachePos, _dataLength - _cachePos);
    this->_fileStream.read(_cache + (_dataLength - _cachePos), _cacheSize - (_dataLength - _cachePos));
    this->_dataLength = this->_fileStream.gcount() + (_dataLength - _cachePos);
    this->_cachePos = 0;
}

char* SortedKMerFile::ReadKmer() {
    if (_accumilated && _cachePos < _dataLength) {
        return GetCurrentKmer();
    }

    char* current = GetCurrentKmer();
    if (current == NULL) {
        return NULL;
    }

    char* next = GetNextKmer();

    while (next != NULL && CheckEquals(current, next)) {

        *(uint32_t*) (next + _kmerStoreLenght - sizeof (uint32_t)) =
                *(uint32_t*) (next + _kmerStoreLenght - sizeof (uint32_t))
                + *(uint32_t*) (current + _kmerStoreLenght - sizeof (uint32_t));

        PopKmer();
        current = GetCurrentKmer();
        next = GetNextKmer();
    }

    this->_accumilated = true;
    return current;
}

void SortedKMerFile::PopKmer() {
    _accumilated = false;
    _cachePos += _kmerStoreLenght;
}

bool SortedKMerFile::CheckEquals(char* lhs, char* rhs) {
    for (int32_t index = 0; index < this->_kmerStoreLenght - sizeof (uint32_t); index += sizeof (uint64_t)) {
        if (*((uint64_t*) (lhs + index)) == *((uint64_t*) (rhs + index))) {
            continue;
        } else {
            return false;
        }
    }
    return true;
}

char* SortedKMerFile::GetCurrentKmer() {
    if (_cachePos + _kmerStoreLenght < _dataLength) {
        return _cache + _cachePos;
    }

    ReadData();

    if (_cachePos < _dataLength) {
        return _cache + _cachePos;
    }
    return NULL;
}

char* SortedKMerFile::GetNextKmer() {
    if (_cachePos + _kmerStoreLenght < _dataLength) {
        return _cache + _cachePos + _kmerStoreLenght;
    }

    ReadData();

    if (_cachePos + _kmerStoreLenght < _dataLength) {
        return _cache + _cachePos + _kmerStoreLenght;
    }
    return NULL;
}