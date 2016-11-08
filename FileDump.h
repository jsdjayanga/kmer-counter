/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FileDump.h
 * Author: jayanga
 *
 * Created on October 21, 2016, 8:26 AM
 */

#pragma once

#include <stdint.h>
#include <string>

#include "FASTQData.h"

using namespace std;

class FileDump {
public:
	FileDump();
    FileDump(string directory);
    FileDump(const FileDump& orig);
    virtual ~FileDump();
    
    void dump(FASTQData* fastqData);
    void dumpToFile(string filename, char* data, int64_t length);
    void dumpKmersToFile(uint32_t id, char* data, int64_t length);
private:
    string _directory;
    int64_t _count;
};


