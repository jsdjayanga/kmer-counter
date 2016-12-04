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
#include <stdint.h>
#include <sstream>
#include <inttypes.h>
#include <thrust/merge.h>
#include "KMerCounter.h"
#include "InputFileHandler.h"
#include "KmerKeyValue.h"

uint32_t __kmer_record_size;

KMerCounter::KMerCounter(Options* options) {
	_options = options;
	_fileDump = new FileDump(_options->getTempFileLocation());
	_gpuStreams = NULL;
	_streamCount = 8;

	uint64_t kmerLength = options->GetKmerLength();
	uint64_t kmerStoreSize = kmerLength / 32;
	if (kmerLength % 32 > 0) {
		kmerStoreSize++;
	}
	kmerStoreSize *= 8;
	kmerStoreSize += 4;

	__kmer_record_size = kmerStoreSize;

	_kMerFileMergeHandler = new KMerFileMergeHandler(_options->getOutputFile(), _options->GetKmerLength(),
			_options->getNoOfMergersAtOnce(), _options->getNoOfMergeThreads());

	_input_done = false;
	_processing_done = false;

	//_countingHashTable = new CountingHashTable<1>(0, 1, 8 * 1000 * 1000 * 1000, 500 * 1000 * 1000);

	_gpuThrustData = PrepareGPUForThrust(200 * 1000 * 1000);
	_host_data_size = 0;
	_host_data_max_size = ((uint64_t)2) * 1024 * 1024 * 1024;

	_host_data_list.clear();

	_temp_file_id = 0;
}

KMerCounter::KMerCounter(const KMerCounter& orig) {
}

KMerCounter::~KMerCounter() {
	delete[] _gpuStreams;
}

void KMerCounter::dispatchWork(GPUStream* gpuStream, FASTQData* fastqData, int64_t lineLength, uint32_t readId) {
	printf("Dispatching work to %i, readid=%i\n", gpuStream, readId);

	//uint64_t outputSize =;
	bool capable = processKMers(gpuStream, _gpuThrustData, fastqData->getData(), _options->GetKmerLength(), fastqData->getSize(),
						lineLength, readId, *_fileDump);
	printf("========================= capable=%i\n", capable);
	while (!capable) {
		printf("======================================waiting for hash table to ready\n");
		gpuStream->_waiting_for_hash_table = true;
		this_thread::sleep_for(chrono::seconds(1));
		capable = processKMers(gpuStream, _gpuThrustData, fastqData->getData(), _options->GetKmerLength(), fastqData->getSize(),
				lineLength, readId, *_fileDump);
		if (capable) {
			gpuStream->_waiting_for_hash_table = false;
		}
	}

//	ostringstream tempFilename;
//	tempFilename << _options->getTempFileLocation() << "/" << readId;
//	_kMerFileMergeHandler->AddFile(tempFilename.str());

//	for (uint64_t i = 0; i < outputSize; i += __kmer_record_size) {
//
//		char* kmer_db = NULL;
//		if (gpuStream->_kmer_db_line_index + __kmer_record_size < gpuStream->_kmer_db_line_length) {
//			kmer_db = gpuStream->_kmer_db.front();
//		} else {
//			gpuStream->_kmer_db.push_front(new char[gpuStream->_kmer_db_line_length]);
//			gpuStream->_kmer_db_line_index = 0;
//			kmer_db = gpuStream->_kmer_db.front();
//		}
//
//		memcpy(kmer_db + gpuStream->_kmer_db_line_index, gpuStream->_h_output + i, __kmer_record_size - sizeof(uint32_t));
//
//		concurrent_hash_map<char*, uint32_t, MyHasher>::accessor acc;
//		if (_con_hashtable.emplace(acc, kmer_db + gpuStream->_kmer_db_line_index,
//				*(uint32_t*)(gpuStream->_h_output + i + __kmer_record_size - sizeof(uint32_t)))) {
//			gpuStream->_kmer_db_line_index += (__kmer_record_size - sizeof(uint32_t));
//		} else {
//			//printf("========================= match found stream:%i\n", gpuStream->_id);
//			acc->second = acc->second + *(uint32_t*)(gpuStream->_h_output + i + __kmer_record_size - sizeof(uint32_t));
//		}
//	}


//	bool capable = false;
//	const uint64_t no_of_records = outputSize / 16;
//	capable = _countingHashTable->Insert((KmerKeyValue<1>*) gpuStream->_d_output, no_of_records);
//	while (!capable) {
//		printf("======================================waiting for hash table to ready\n");
//		gpuStream->_waiting_for_hash_table = true;
//		this_thread::sleep_for(chrono::milliseconds(50));
//		capable = _countingHashTable->Insert((KmerKeyValue<1>*) gpuStream->_d_output, no_of_records);
//		if (capable) {
//			gpuStream->_waiting_for_hash_table = false;
//		}
//	}

//	printf("========================= outsize: %"PRIu64", capable=%i\n", capable);

	_rec_mtx.lock();
	_vacantStreams.insert(gpuStream);
	_rec_mtx.unlock();

	delete fastqData;
}

void KMerCounter::DumpResults() {
	while (!_processing_done) {
		this_thread::sleep_for(chrono::seconds(1));
		//_countingHashTable->UpdateCounters();
//		printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++%"PRIu64"\n", _countingHashTable->GetSuccessCount());
		printf("++++++++++++++++++++++++++++++++++++++++++++++++++DumpResults called+++++++\n");

		int32_t count = 0;
		for (int32_t i = 0; i < _streamCount; i++) {
			if (_gpuStreams[i]->_waiting_for_hash_table) {
				count++;
			}
		}

		if (count >= _streamCount || _input_done) {
			printf("======================================Need to to temp sort and count\n");

			ThrustProcessedResut* data = processThrust(_gpuThrustData, 32);
//			uint64_t count = 0;
//			memcpy(&count, data, sizeof(uint64_t));
			_host_data_size += data->_count * sizeof(KMer32) + data->_count * sizeof(uint32_t);
			_host_data_list.push_back(data);

//			std::ofstream output_file("/tmp/1/testoutput.log");
//			char* data = new char[_gpuThrustData->_length];
//			CUDA_CHECK_RETURN(cudaMemcpy(data, ((char*)_gpuThrustData->_d_data), _gpuThrustData->_length, cudaMemcpyDeviceToHost));
//			for (int64_t j = 0; j < _gpuThrustData->_length; j += 8) {
//					output_file.write((data + j), 8);
//			}
//
//			delete[] data;
//			output_file.close();

//			printf("======================================Dumping files\n");
//			_countingHashTable->TempDump();
//			_countingHashTable->Init();
//			printf("======================================Dumping files completed");
		}

		if (!(_host_data_size < _host_data_max_size)) {
			list<ThrustProcessedResut*> temp_list;
			temp_list.splice(temp_list.begin(), _host_data_list);
			MergeData(temp_list, true);
		}

	}

	ThrustProcessedResut* data = processThrust(_gpuThrustData, 32);
	_host_data_size += data->_count * sizeof(KMer32) + data->_count * sizeof(uint32_t);
	_host_data_list.push_back(data);

	list<ThrustProcessedResut*> temp_list;
	temp_list.splice(temp_list.begin(), _host_data_list);

	if (_temp_file_id > 0) {
		MergeData(temp_list, true);

		// TODO : need to merge files;
	} else {
		MergeData(temp_list, false);
		ofstream output_file(_options->getOutputFile());
		ThrustProcessedResut* data = _host_data_list.front();
		_host_data_list.pop_front();
		for (uint64_t i = 0; i < data->_count; i++) {
			//std::cout << " " << it->first << ":" << it->second;
			if (*((uint32_t*)data->_values + i) != 0) {
				output_file.write((char*)(((KMer32*)data->_keys) + i), sizeof(KMer32));
				output_file.write((char*)(((uint32_t*)data->_values) + i), sizeof(uint32_t));
			}
		}
		output_file.close();
		delete data;
	}


//	ofstream output_file(_options->getOutputFile());
//	for (auto it = _con_hashtable.begin(); it != _con_hashtable.end(); ++it) {
////	for (auto it = _con_uo_hashtable.begin(); it != _con_uo_hashtable.end(); ++it) {
//		//std::cout << " " << it->first << ":" << it->second;
//		output_file.write(it->first, 8);
//		output_file.write((char*)&(it->second), 4);
//	}
//	output_file.close();

//	_countingHashTable->Dump();
}

void KMerCounter::MergeData(list<ThrustProcessedResut*> data_list, bool write_to_file) {
	while(data_list.size() > 1) {
		ThrustProcessedResut* data1 = data_list.front();
		data_list.pop_front();

		ThrustProcessedResut* data2 = data_list.front();
		data_list.pop_front();

//		thrust::host_vector key1((KMer32*)(data1 + sizeof(uint64_t)), (KMer32*)(data1 + sizeof(uint64_t)) + count1);
//		thrust::host_vector value1((uint32_t*)(data1 + sizeof(uint64_t) + count * sozeof(KMer32)), (uint32_t*)(data1 + sizeof(uint64_t) + count * sozeof(KMer32)) + count1);
//
//		thrust::host_vector key2((KMer32*)(data2 + sizeof(uint64_t)), (KMer32*)(data2 + sizeof(uint64_t)) + count2);
//		thrust::host_vector value2((uint32_t*)(data2 + sizeof(uint64_t) + count * sozeof(KMer32)), (uint32_t*)(data2 + sizeof(uint64_t) + count * sozeof(KMer32)) + count1);

		char* temp1 = new char[(data1->_count + data2->_count) * sizeof(KMer32)];
		char* temp2 = new char[(data1->_count + data2->_count) * sizeof(uint32_t)];
		thrust::merge_by_key((KMer32*)(data1->_keys), (KMer32*)(data1->_keys) + data1->_count,
				             (KMer32*)(data2->_keys), (KMer32*)(data2->_keys) + data2->_count,
				             (uint32_t*)(data1->_values),
				             (uint32_t*)(data2->_values),
				             (KMer32*)temp1,
				             (uint32_t*)temp2,
				             lessthan_key());

		thrust::host_vector<KMer32> keys_out(data1->_count + data2->_count);
		thrust::host_vector<int32_t> lengths(data1->_count + data2->_count);
		int rsize = thrust::reduce_by_key((KMer32*)temp1,
				((KMer32*)temp1) + data1->_count + data2->_count,
				(uint32_t*)temp2,
					keys_out.begin(),
					lengths.begin(),
					equal_key()).first - keys_out.begin();

		uint64_t count = keys_out.size();
		char* temp_keys = new char[count * sizeof(KMer32)];
		char* temp_values = new char[count * sizeof(uint32_t)];

		thrust::copy(keys_out.begin(), keys_out.end(), (KMer32*)temp_keys);
		thrust::copy(lengths.begin(), lengths.end(), (uint32_t*)temp_values);
		ThrustProcessedResut* tps = new ThrustProcessedResut(count, temp_keys, temp_values);

		delete[] temp1;
		delete[] temp2;
		delete data1;
		delete data2;

		data_list.push_back(tps);
	}

	if (write_to_file == true) {
		//TODO: write to temp file
		string temp_file_name = _options->getTempFileLocation() + "/" + to_string(_temp_file_id++);
		ofstream output_file(temp_file_name.c_str());
		ThrustProcessedResut* data = data_list.front();
		data_list.pop_front();
		for (uint64_t i = 0; i < data->_count; i++) {
			//std::cout << " " << it->first << ":" << it->second;
			if (*((uint32_t*)data->_values + i) != 0) {
				output_file.write((char*)(((KMer32*)data->_keys) + i), sizeof(KMer32));
				output_file.write((char*)(((uint32_t*)data->_values) + i), sizeof(uint32_t));
			}
		}
		output_file.close();
		delete data;

	} else {
		_host_data_list.splice(_host_data_list.begin(), data_list);
	}
}

void KMerCounter::Start() {
	uint32_t readId = 0;
	InputFileHandler* inputFileHandler = new InputFileHandler(_options->GetInputFileDirectory());
	//_kMerFileMergeHandler->Start();
	thread finalizer(&KMerCounter::DumpResults, this);
	int64_t chunkSize = GetChunkSize(inputFileHandler->getLineLength(), _options->GetKmerLength(),
			_options->GetGpuMemoryLimit());
	FASTQData* fastqData = inputFileHandler->read(chunkSize);

	_gpuStreams = PrepareGPU(_streamCount, fastqData->getSize(), inputFileHandler->getLineLength(), _options->GetKmerLength());
	for (int i = 0; i < _streamCount; i++) {
		_vacantStreams.insert(_gpuStreams[i]);
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

	_input_done = true;

	for (std::list<thread>::iterator it=_workerThreads.begin(); it != _workerThreads.end(); ++it) {
	    (*it).join();;
	}
	_processing_done = true;
	finalizer.join();

	FreeGPU(_gpuStreams, _streamCount);

	for (int i = 0; i < _streamCount; i++) {
		GPUStream* stream = _gpuStreams[i];
		for (list<char*>::iterator it=stream->_kmer_db.begin(); it != stream->_kmer_db.end(); ++it) {
			delete[] (*it);
		}
		delete stream;
	}
	delete inputFileHandler;

	// Count KMers with Merged Files
//	_kMerFileMergeHandler->InputComplete();
//	_kMerFileMergeHandler->Join();

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
