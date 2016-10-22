/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FASTQFileReader.h
 * Author: jayanga
 *
 * Created on October 19, 2016, 5:45 PM
 */

#pragma once

#include <string>
#include <fstream>

#include "FASTQData.h"

using namespace std;

class FASTQFileReader {
public:
    FASTQFileReader(string filename, int64_t fileSize);
    FASTQFileReader(const FASTQFileReader& orig);
    virtual ~FASTQFileReader();
    
    FASTQData* readData(int64_t chunkSize);
    bool isComplete();
    int64_t getLineLength();
private:
    ifstream _fileStream;
    string _filename;
    int64_t _fileSize;
    int64_t _lineLength;
    int64_t _chunkSize;
};


