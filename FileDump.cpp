/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FileDump.cpp
 * Author: jayanga
 * 
 * Created on October 21, 2016, 8:26 AM
 */

#include <fstream>
#include <algorithm>
#include <string>
#include <sstream>

#include "FileDump.h"

FileDump::FileDump(string directory) {
    this->_directory = directory;
    this->_count = 0;
}

FileDump::FileDump(const FileDump& orig) {
}

FileDump::~FileDump() {
}

void FileDump::dump(FASTQData* fastqData) {
    ofstream output;
    ostringstream file;
    file << _directory << "/" << _count++;
    output.open(file.str().c_str());
    output.write(fastqData->getData(), fastqData->getSize());
    output.close();
}
