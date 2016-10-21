/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   KMerCounter.cpp
 * Author: jayanga
 * 
 * Created on October 19, 2016, 8:03 PM
 */

#include "KMerCounter.h"
#include "InputFileHandler.h"

KMerCounter::KMerCounter(Options* options) {
    _options = options;
    _fileDump = new FileDump("/tmp/1");
}

KMerCounter::KMerCounter(const KMerCounter& orig) {
}

KMerCounter::~KMerCounter() {
}

void KMerCounter::Start() {
    InputFileHandler* inputFileHandler = new InputFileHandler(_options->GetInputFileDirectory(), _options->GetChunkSize());
    FASTQData* fastqData = inputFileHandler->read();
    while (fastqData != NULL) {
        // TODO : Pump data to GPU
        //cout << "====" << fastqData->getLineLength() << endl;
        _fileDump->dump(fastqData);
        
        fastqData = inputFileHandler->read();
    }

    
//    list<InputFileDetails*>& list = inputFileHandler->getFileList();
    
//    for (std::list<InputFileDetails*>::iterator it = list.begin(); it != list.end(); it++) {
//        cout << (*it)->GetFilename() << endl;
//    }
}