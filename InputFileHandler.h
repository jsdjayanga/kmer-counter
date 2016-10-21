/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   InputFileHandler.h
 * Author: jayanga
 *
 * Created on October 19, 2016, 8:20 PM
 */

#pragma once

#include <string>
#include <list>
#include "FASTQData.h"
#include "FASTQFileReader.h"

using namespace std;

class InputFileHandler {
public:
    InputFileHandler(string fileDirectory, int64_t chunkSize);
    InputFileHandler(const InputFileHandler& orig);
    virtual ~InputFileHandler();
    
    //list<InputFileDetails*>& getFileList();
    FASTQData* read();
private:
    string _fileDirectory;
    int64_t _chunkSize;
    list<FASTQFileReader*> _fileReaders;
};


