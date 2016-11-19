/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   KMerFileMergeHandler.h
 * Author: jayanga
 *
 * Created on November 17, 2016, 6:08 PM
 */

#pragma once

#include <thread>
#include <mutex>
#include <list>
#include <string>
#include <vector>
#include <set>
#include <stdint.h>

#include "KMerFileMerger.h"

using namespace std;

class KMerFileMergeHandler {
public:
    KMerFileMergeHandler(string outputFilename, uint64_t kmerLength, uint32_t noOfMergersAtOnce, uint32_t maxThreads);
    KMerFileMergeHandler(const KMerFileMergeHandler& orig);
    virtual ~KMerFileMergeHandler();
    void Start();
    void Join();
    void AddFile(string filename);
    void InputComplete();
private:
    string _outputFilename;
    thread _mergerThread;
    list<thread> _workerThreads;
    recursive_mutex _rec_mtx;
    list<string> _files;
    set<KMerFileMerger*> _kmerMergers;
    bool _inputComplete;
    bool _completed;
    uint32_t _maxThreads;
    
    uint32_t _noOfMergersAtOnce;
    uint64_t _kmerLength;

    void Run();
    void PerformMerge(KMerFileMerger* kMerFileMerger, string outputFilename);
};


