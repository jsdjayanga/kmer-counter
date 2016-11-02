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
#include "GPUHandler.h"

KMerCounter::KMerCounter(Options* options) {
    _options = options;
    _fileDump = new FileDump("/tmp/1");
}

KMerCounter::KMerCounter(const KMerCounter& orig) {
}

KMerCounter::~KMerCounter() {
}

void KMerCounter::Start() {
    InputFileHandler* inputFileHandler = new InputFileHandler(_options->GetInputFileDirectory());
    int64_t chunkSize = GetChunkSize(inputFileHandler->getLineLength(), _options->GetKmerLength(), _options->GetGpuMemoryLimit());
    FASTQData* fastqData = inputFileHandler->read(chunkSize);
    while (fastqData != NULL) {
        // TODO : Pump data to GPU
        //cout << "====" << fastqData->getLineLength() << endl;
        _fileDump->dump(fastqData);
        
        processKMers(fastqData->getData(), _options->GetKmerLength(), fastqData->getSize(), inputFileHandler->getLineLength());
        
        chunkSize = GetChunkSize(inputFileHandler->getLineLength(), _options->GetKmerLength(), _options->GetGpuMemoryLimit());
        fastqData = inputFileHandler->read(chunkSize);
    }

    
//    list<InputFileDetails*>& list = inputFileHandler->getFileList();
    
//    for (std::list<InputFileDetails*>::iterator it = list.begin(); it != list.end(); it++) {
//        cout << (*it)->GetFilename() << endl;
//    }
}

int64_t KMerCounter::GetChunkSize(int64_t lineLength, int64_t kmerLength, int64_t gpuMemoryLimit) {
    int64_t bytesNeededForKmer = kmerLength / 4;
    if (kmerLength % 4 != 0) {
        bytesNeededForKmer++;
    }
    
    int64_t byteRepresentationForKmer = bytesNeededForKmer / 8;
    if (bytesNeededForKmer % 8 != 0) {
        byteRepresentationForKmer++;
    }
    byteRepresentationForKmer++; // for the k-mer count
    byteRepresentationForKmer *= 8;
    
    int64_t numberOfKmers = (lineLength - kmerLength + 1);
    int64_t memoryForAllKmers = byteRepresentationForKmer * numberOfKmers;
    
    int64_t possibleRecordCount = (gpuMemoryLimit - lineLength) / (memoryForAllKmers - 1);

    return lineLength * possibleRecordCount;
}
