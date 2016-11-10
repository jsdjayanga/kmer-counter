/*
 * KMerFileReader.cpp
 *
 *  Created on: Nov 10, 2016
 *      Author: jayanga
 */

#include <string.h>
#include "KMerFileReader.h"

KMerFileReader::KMerFileReader(string filename, uint64_t fileSize, uint64_t recordLength, uint64_t recordCount) {
	// TODO Auto-generated constructor stub
	this->_filename = filename;
	this->_fileSize = fileSize;
	this->_recordLength = recordLength;
	this->_recordCount = recordCount;

	this->_cacheSize = recordLength * recordCount;
	this->_cache = new char[this->_cacheSize];
	memset(this->_cache, 0, this->_cacheSize);
	this->_cachePos = 0;

	_fileStream.open(_filename.c_str());
	if (_fileStream.is_open()) {
		_fileStream.seekg(0);
	}

	_fileStream.read(_cache, _cacheSize);
	_availableLength = _fileStream.gcount();
	for (uint64_t index = 0; index < _availableLength; index += _recordLength) {
		uint32_t* count = (uint32_t*) (_cache + index + _recordLength - sizeof(uint32_t));
		if (*count > 0) {
			memmove(_cache, _cache + index, _cacheSize - index);
			_fileStream.read(_cache + (_cacheSize - index), index);
			_availableLength = _availableLength - index + _fileStream.gcount();
			break;
		}
		if (index + _recordLength >= _availableLength) {
			_fileStream.read(_cache, _cacheSize);
			_availableLength = _fileStream.gcount();
			index = 0;
		}
	}
}

KMerFileReader::~KMerFileReader() {
	// TODO Auto-generated destructor stub
}

char* KMerFileReader::peekKmer() {
    if (_cachePos < _availableLength) {
        return _cache + _cachePos;
    } else {
        readData(0, _cacheSize);
        if (_cachePos < _availableLength) {
            return _cache + _cachePos;
        }
    }
    return NULL;
}

char* KMerFileReader::peekNextKmer() {
    if (_cachePos < _availableLength - _recordLength) {
        return _cache + _cachePos + _recordLength;
    } else {
        readData(0, _cacheSize);
        if (_cachePos < _availableLength - _recordLength) {
            return _cache + _cachePos;
        }
    }
    return NULL;
}

char* KMerFileReader::popKmer() {
    if (_cachePos < _availableLength) {
        uint64_t currentPos = _cachePos;
        _cachePos += _recordLength;
        return _cache + currentPos;
    } else {
        memmove(_cache, _cache + _cachePos, _cacheSize - _cachePos);
        readData(_recordLength, _cacheSize - _recordLength);
        return _cache;
    }
}

void KMerFileReader::readData(uint64_t cacheStartPos, uint64_t size) {
    _fileStream.read(_cache + cacheStartPos, size);
    if (_fileStream.gcount() > 0) {
        _availableLength = _fileStream.gcount() + _cacheSize - _cachePos;
        _cachePos = 0;
    }
}
