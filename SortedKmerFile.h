/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   SortedKmerFile.h
 * Author: jayanga
 *
 * Created on December 8, 2016, 6:37 PM
 */

#ifndef SORTEDKMERFILE_H
#define SORTEDKMERFILE_H

#include <string>
#include <fstream>
#include <stdint.h>

class SortedKmerFile {
public:
    SortedKmerFile(std::string filename, uint32_t kmerLength);
    SortedKmerFile(const SortedKmerFile& orig);
    virtual ~SortedKmerFile();
    char* Peek();
    void UpdatePeek(uint32_t count);
    void Pop();
private:
    std::string _filename;
    std::ifstream _file_stream;
    char* _buffer;
    uint32_t _buffer_size;
    uint32_t _buffer_index;
    uint32_t _buffer_valid_index;
    uint32_t _kmerLength;
    uint32_t _kmerStoreLength;
    
    void Read();
};

#endif /* SORTEDKMERFILE_H */

