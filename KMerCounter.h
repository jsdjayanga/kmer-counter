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

#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_unordered_map.h>
#include <iostream>
#include <thread>
#include <mutex>
#include "Options.h"
#include "FileDump.h"
#include "KMerFileMerger.h"
#include "KMerFileMergeHandler.h"
#include "GPUHandler.h"
#include <unordered_map>


using namespace std;
using namespace tbb;

struct eqstr
{
	bool operator()(const char* s1, const char* s2) const
	{
		return (s1 == s2) || (s1 && s2 && *(uint64_t*)s1 == *(uint64_t*)s2);
	}
};

struct MyHasher {
	static size_t hash(const char* value)                  { return *(uint64_t*)value; }
	static bool   equal(const char* s1, const char* s2) { return (s1 == s2) || (s1 && s2 && *(uint64_t*)s1 == *(uint64_t*)s2); }
};

class KMerCounter {
public:
    KMerCounter(Options* options);
    KMerCounter(const KMerCounter& orig);
    virtual ~KMerCounter();
    void Start();
private:
    Options* _options;
    FileDump* _fileDump;
    KMerFileMergeHandler* _kMerFileMergeHandler;
    list<thread> _workerThreads;
    recursive_mutex _rec_mtx;
    set<GPUStream*> _vacantStreams;
    
    int64_t GetChunkSize(int64_t lineLength, int64_t kmerLength, int64_t gpuMemoryLimit);
    void dispatchWork(GPUStream* gpuStream, FASTQData* fastqData, int64_t lineLength, uint32_t readId);
    void DumpResults();


    bool _processing_done;

    recursive_mutex _rec_mtxs[8];

    concurrent_hash_map<char*, uint32_t, MyHasher> _con_hashtable;
    //concurrent_unordered_map<char*, uint32_t, MyHasher, eqstr> _con_uo_hashtable;
};


