/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   KMerCounter.h
 * Author: jayanga
 *
 * Created on October 19, 2016, 8:03 PM
 */

#pragma once

#include <iostream>
#include "Options.h"
#include "FileDump.h"
#include "KMerFileMerger.h"

using namespace std;

class KMerCounter {
public:
    KMerCounter(Options* options);
    KMerCounter(const KMerCounter& orig);
    virtual ~KMerCounter();
    void Start();
private:
    Options* _options;
    FileDump* _fileDump;
    KMerFileMerger* _kMerFileMerger;
    
    int64_t GetChunkSize(int64_t lineLength, int64_t kmerLength, int64_t gpuMemoryLimit);
};


