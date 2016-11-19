/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   KMerFileMergeHandler.cpp
 * Author: jayanga
 * 
 * Created on November 17, 2016, 6:08 PM
 */

#include <unistd.h>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <set>
#include <cmath>
#include "KMerFileMergeHandler.h"
#include "KMerFileMerger.h"

KMerFileMergeHandler::KMerFileMergeHandler(string outputFilename, uint64_t kmerLength, uint32_t noOfMergersAtOnce, uint32_t maxThreads) {
    this->_outputFilename = outputFilename;
    this->_kmerLength = kmerLength;
    this->_maxThreads = maxThreads;
    this->_workerThreads.clear();
    this->_kmerMergers.clear();
    this->_files.clear();
    this->_noOfMergersAtOnce = noOfMergersAtOnce;
    _inputComplete = false;
    _completed = false;
}

KMerFileMergeHandler::KMerFileMergeHandler(const KMerFileMergeHandler& orig) {
}

KMerFileMergeHandler::~KMerFileMergeHandler() {
}

void KMerFileMergeHandler::Start() {
    _mergerThread = thread(&KMerFileMergeHandler::Run, this);
}

void KMerFileMergeHandler::Join() {
    _mergerThread.join();
}

void KMerFileMergeHandler::Run() {
    uint32_t index = 0;

    while (!_completed) {
        
        if (_kmerMergers.size() >= _maxThreads || _noOfMergersAtOnce > _files.size()) {
            this_thread::sleep_for(chrono::seconds(1));
            continue;
        }
        
        _rec_mtx.lock();

        if (_kmerMergers.size() < _maxThreads && _noOfMergersAtOnce <= _files.size()) {
            list<string> filesToMerge;
            int32_t i = 0;
            for (; i < _noOfMergersAtOnce; i++) {
                filesToMerge.push_back(_files.front());
                _files.pop_front();
            }

            if (filesToMerge.empty()) {
                _rec_mtx.unlock();
                continue;
            }

            string outputFileName = filesToMerge.front() + "_" + to_string(index++);
            KMerFileMerger* kMerFileMerger = new KMerFileMerger(filesToMerge, outputFileName, _kmerLength);

            _kmerMergers.insert(kMerFileMerger);
            _workerThreads.push_back(thread(&KMerFileMergeHandler::PerformMerge, this, kMerFileMerger, outputFileName));
        }

        if (_inputComplete && _kmerMergers.size() + _files.size() / _noOfMergersAtOnce <= 1) {
            _completed = true;
        }

        _rec_mtx.unlock();

    }

    for (std::list<thread>::iterator it = _workerThreads.begin(); it != _workerThreads.end(); ++it) {
        (*it).join();
    }

    list<string> filesToMerge;
    for (std::list<string>::iterator it = _files.begin(); it != _files.end(); ++it) {
        filesToMerge.push_back(*it);
    }
    KMerFileMerger* kMerFileMerger = new KMerFileMerger(filesToMerge, _outputFilename, _kmerLength);
    kMerFileMerger->Merge();
    delete kMerFileMerger;
}

void KMerFileMergeHandler::AddFile(string filename) {
    _rec_mtx.lock();
    _files.push_back(filename);
    _rec_mtx.unlock();
}

void KMerFileMergeHandler::PerformMerge(KMerFileMerger* kMerFileMerger, string outputFilename) {
    kMerFileMerger->Merge();

    _rec_mtx.lock();
    _files.push_back(outputFilename);
    _kmerMergers.erase(kMerFileMerger);
    _rec_mtx.unlock();
    
    delete kMerFileMerger;
}

void KMerFileMergeHandler::InputComplete() {
    _rec_mtx.lock();
    _inputComplete = true;
    _rec_mtx.unlock();
}
