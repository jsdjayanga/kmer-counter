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

#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include "KMerCounter.h"
#include "InputFileHandler.h"


KMerCounter::KMerCounter(Options* options) {
	_options = options;
	_fileDump = new FileDump(_options->getTempFileLocation());

	uint64_t kmerLength = options->GetKmerLength();
	uint64_t kmerStoreSize = kmerLength / 32;
	if (kmerLength % 32 > 0) {
		kmerStoreSize++;
	}
	kmerStoreSize *= 8;
	kmerStoreSize += 4;

	_kMerFileMergeHandler = new KMerFileMergeHandler(_options->getOutputFile(), _options->GetKmerLength(),
			_options->getNoOfMergersAtOnce(), _options->getNoOfMergeThreads());
}

KMerCounter::KMerCounter(const KMerCounter& orig) {
}

KMerCounter::~KMerCounter() {
}

void KMerCounter::dispatchWork(GPUStream* gpuStream, FASTQData* fastqData, int64_t lineLength, uint32_t readId) {
	printf("Dispatching work to %i, readid=%i", gpuStream, readId);

	processKMers(gpuStream, fastqData->getData(), _options->GetKmerLength(), fastqData->getSize(),
						lineLength, readId, *_fileDump);

	ostringstream tempFilename;
	tempFilename << _options->getTempFileLocation() << "/" << readId;
	_kMerFileMergeHandler->AddFile(tempFilename.str());

	_rec_mtx.lock();
	_vacantStreams.insert(gpuStream);
	_rec_mtx.unlock();

	delete fastqData;
}

void KMerCounter::Start() {
	uint32_t readId = 0;
	InputFileHandler* inputFileHandler = new InputFileHandler(_options->GetInputFileDirectory());
	_kMerFileMergeHandler->Start();
	int64_t chunkSize = GetChunkSize(inputFileHandler->getLineLength(), _options->GetKmerLength(),
			_options->GetGpuMemoryLimit());
	FASTQData* fastqData = inputFileHandler->read(chunkSize);

	uint32_t streamCount = 8;
	GPUStream** gpuStreams = PrepareGPU(streamCount, fastqData->getSize(), inputFileHandler->getLineLength(), _options->GetKmerLength());
	for (int i = 0; i < streamCount; i++) {
		_vacantStreams.insert(gpuStreams[i]);
	}

	while (fastqData != NULL) {

		readId++;
		// TODO : Pump data to GPU
		//cout << "====" << fastqData->getLineLength() << endl;
		//_fileDump->dump(fastqData);

		if (fastqData->getSize() > 0 && fastqData->getSize() >= inputFileHandler->getLineLength()) {
			while (_vacantStreams.empty()) {
				this_thread::sleep_for(chrono::milliseconds(50));
			}

			_rec_mtx.lock();
			GPUStream* gpuStream = *(_vacantStreams.begin());
			_vacantStreams.erase(_vacantStreams.begin());
			_workerThreads.push_back(thread(&KMerCounter::dispatchWork, this, gpuStream, fastqData, inputFileHandler->getLineLength(), readId));
			_rec_mtx.unlock();
		}

		fastqData = inputFileHandler->read(chunkSize);
	}

	for (std::list<thread>::iterator it=_workerThreads.begin(); it != _workerThreads.end(); ++it) {
	    (*it).join();;
	}

	// Count KMers with Merged Files
	_kMerFileMergeHandler->InputComplete();
	_kMerFileMergeHandler->Join();

//	string tempLocation = _options->getTempFileLocation();
//	list<string> tempFiles;
//	DIR *dir;
//	struct dirent *ent;
//	if ((dir = opendir(tempLocation.c_str())) != NULL) {
//		while ((ent = readdir(dir)) != NULL) {
//			if (strncmp(".", ent->d_name, 1) != 0 && strncmp("..", ent->d_name, 2) != 0) {
//				ifstream temp;
//				string filename = tempLocation + "/" + ent->d_name;
//				tempFiles.push_back(filename);
//			}
//		}
//		closedir(dir);
//	} else {
//		cout << "Couldn't open directory : " << tempLocation << endl;
//	}
//	KMerFileMerger* kmerMerger = new KMerFileMerger(tempFiles, _options->getOutputFile(), _options->GetKmerLength());
//	kmerMerger->Merge();

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
