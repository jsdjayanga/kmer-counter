/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   SortedKmerFile.cpp
 * Author: jayanga
 * 
 * Created on December 8, 2016, 6:37 PM
 */

#include <string.h>
#include <iostream>
#include "SortedKmerFile.h"

using namespace std;

SortedKmerFile::SortedKmerFile(std::string filename, uint32_t kmerLength) {
    this->_filename = filename;
    this->_file_stream.open(this->_filename.c_str());
    this->_kmerLength = kmerLength;
    this->_kmerStoreLength = (((this->_kmerLength + 31) / 32) * sizeof (uint64_t)) + sizeof (uint32_t);

    this->_buffer_size = 12 * 1024 * 1024;
    this->_buffer_index = 0;
    this->_buffer_valid_index = 0;

    this->_buffer = new char[_buffer_size];
    
    if (this->_file_stream.is_open()) {
        Read();
    }
}

SortedKmerFile::SortedKmerFile(const SortedKmerFile& orig) {
}

SortedKmerFile::~SortedKmerFile() {
}

void SortedKmerFile::Read() {
//    cout << "=========Read" << endl;
    
    
    if (this->_buffer_index > 0) {
        memmove(this->_buffer, this->_buffer + this->_buffer_index, this->_buffer_valid_index - this->_buffer_index);
        this->_file_stream.read(this->_buffer + (this->_buffer_valid_index - this->_buffer_index), this->_buffer_size - (this->_buffer_valid_index - this->_buffer_index));
        this->_buffer_valid_index = _file_stream.gcount() + (this->_buffer_valid_index - this->_buffer_index);
    } else {
        this->_file_stream.read(this->_buffer, this->_buffer_size);
        this->_buffer_valid_index = _file_stream.gcount();
    }
    this->_buffer_index = 0;
}

char* SortedKmerFile::Peek() {
//    cout << "=========Peek, " << this->_buffer_index << "|" << this->_buffer_valid_index << "|" << this->_buffer_size  << endl;
    
    if (this->_buffer_index < this->_buffer_valid_index) {
        return _buffer + this->_buffer_index;
    } else {
        Read();
        if (this->_buffer_index < this->_buffer_valid_index) {
            return _buffer + this->_buffer_index;
        } else {
            return NULL;
        }
    }
}

void SortedKmerFile::UpdatePeek(uint32_t count) {
    *(uint32_t*) (_buffer + this->_buffer_index + this->_kmerStoreLength - sizeof (uint32_t)) += count;
}

void SortedKmerFile::Pop() {
//    cout << "=========Pop, " << this->_buffer_index << "|" << this->_buffer_valid_index << "|" << this->_buffer_size  << endl;
    this->_buffer_index += this->_kmerStoreLength;
}