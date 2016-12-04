/*
 * CountingHashTable.h
 *
 *  Created on: Nov 29, 2016
 *      Author: jayanga
 */

#ifndef COUNTINGHASHTABLEBASE_H_
#define COUNTINGHASHTABLEBASE_H_

#include <vector>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
//#include <stdexcept>
#include <inttypes.h>
#include "KmerKeyValue.h"

/**
 * This macro checks return value of the CUDA runtime call and exits
 * the application if the call failed.
 */
#define CUDA_CHECK_RETURN(value) {											\
	cudaError_t _m_cudaStat = value;										\
	if (_m_cudaStat != cudaSuccess) {										\
		fprintf(stderr, "Error %s at line %d in file %s\n",					\
				cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);		\
		exit(1);															\
	} }



//template<class KeyType>
//constexpr uint32_t key_size() {
//	return sizeof(KeyType) / sizeof(uint64_t);
//}

template<uint32_t key_size>
extern void InsertToHashTable(KmerKeyValue<key_size>* input, uint32_t no_of_keys_per_stream, cudaStream_t stream,
		KmerKeyValue<key_size>* kmer_db, uint64_t kmer_db_max_record_count, uint64_t* cuda_counters);


template<uint32_t key_size>
class CountingHashTableBase {
private:
	char* _d_data_buffer;
	uint64_t _kmer_db_size;
	int64_t _kmer_db_record_count;
	int64_t _kmer_db_max_record_count;
	KmerKeyValue<key_size>* _d_kmer_db;

	char* _h_failed_data_buffer;
	char* _d_failed_data_buffer;
	uint64_t _failed_kmer_size;
	int64_t _failed_kmer_entries_count;
	int64_t _failed_kmer_max_entry_count;
	KmerKeyValue<key_size>* _d_failed_kmer_entires;

	// Environment information
	uint32_t _device_id;
	uint32_t _no_of_streams;
	cudaStream_t* _cuda_streams;

	uint64_t* _h_cuda_counters;
	uint64_t* _d_cuda_counters;
	uint64_t _cuda_counters_per_stream;

public:
	CountingHashTableBase(const uint32_t device_id, const uint32_t no_of_streams, const uint32_t kmer_db_size, const uint32_t kmer_failed_db_size) :
			_device_id(device_id), _d_data_buffer(NULL), _d_kmer_db(NULL), _kmer_db_size(kmer_db_size), _no_of_streams(no_of_streams), _kmer_db_record_count(
					0), _kmer_db_max_record_count(0), _failed_kmer_size(kmer_failed_db_size), _failed_kmer_entries_count(0), _failed_kmer_max_entry_count(
					0) {

		// set device id
		CUDA_CHECK_RETURN(cudaSetDevice(_device_id));

		// determining usable free memory
//		uint64_t free_memory, total_memory;
//		cudaMemGetInfo(&free_memory, &total_memory);
//		_kmer_db_size = 0.75 * free_memory;
//		_failed_kmer_size = 10 * 1024 * 1024;

		_cuda_counters_per_stream = 2;

		// allocating memory for _kmer_db
		CUDA_CHECK_RETURN(cudaMalloc((void ** ) &_d_data_buffer, _kmer_db_size));
		CUDA_CHECK_RETURN(cudaMalloc((void ** ) &_d_failed_data_buffer, _failed_kmer_size));

		_h_cuda_counters = new uint64_t[_cuda_counters_per_stream * _no_of_streams];
		CUDA_CHECK_RETURN(cudaMalloc((void ** ) &_d_cuda_counters, _cuda_counters_per_stream * _no_of_streams * sizeof(uint64_t)));

		_d_kmer_db = (KmerKeyValue<key_size>*)_d_data_buffer;
		_d_failed_kmer_entires = (KmerKeyValue<key_size>*)_d_failed_data_buffer;

//		_h_failed_data_buffer = new char[_failed_kmer_size];

		// create streams
		_cuda_streams = new cudaStream_t[_no_of_streams];
		for (uint32_t i = 0; i < _no_of_streams; i++) {
			CUDA_CHECK_RETURN(cudaStreamCreate(&_cuda_streams[i]));
		}

		// other
		_kmer_db_max_record_count = _kmer_db_size / sizeof(KmerKeyValue<key_size> );
		_failed_kmer_max_entry_count = _failed_kmer_size / sizeof(KmerKeyValue<key_size>);

		Init();
	}

	virtual ~CountingHashTableBase() {
		// set device id
		CUDA_CHECK_RETURN(cudaSetDevice(_device_id));

		// free allocated memory
		CUDA_CHECK_RETURN(cudaFree(_d_data_buffer));
		CUDA_CHECK_RETURN(cudaFree(_d_failed_data_buffer));

		// destroy streams
		for (uint32_t i = 0; i < _no_of_streams; i++) {
			CUDA_CHECK_RETURN(cudaStreamDestroy(_cuda_streams[i]));
		}
	}

	void Init() {
		// set device id
		CUDA_CHECK_RETURN(cudaSetDevice(_device_id));
		CUDA_CHECK_RETURN(cudaDeviceSynchronize());

//		CUDA_CHECK_RETURN(cudaMemcpy(_h_failed_data_buffer, _d_failed_data_buffer, _failed_kmer_size, cudaMemcpyDeviceToHost));
//		UpdateCounters();

		// reset allocated memory
		CUDA_CHECK_RETURN(cudaMemset(_d_data_buffer, 0, _kmer_db_size));
//		CUDA_CHECK_RETURN(cudaMemset(_d_failed_data_buffer, 0, _failed_kmer_size));
		CUDA_CHECK_RETURN(cudaMemset(_d_cuda_counters, 0, _cuda_counters_per_stream * _no_of_streams * sizeof(uint64_t)));

		UpdateCounters();
		if (_failed_kmer_entries_count > 0) {
			Insert((KmerKeyValue<key_size>*)_d_failed_data_buffer, _failed_kmer_entries_count);
		}
		CUDA_CHECK_RETURN(cudaMemset(_d_failed_data_buffer, 0, _failed_kmer_size));

		_kmer_db_record_count = 0;
		_failed_kmer_entries_count = 0;
	}

	bool Insert(KmerKeyValue<key_size>* d_input, const uint64_t no_of_keys) {
		// set device id
		CUDA_CHECK_RETURN(cudaSetDevice(_device_id));

		// TODO : find the proper ratio for this
		if (no_of_keys * 1.25 > _kmer_db_max_record_count - _kmer_db_record_count
				|| _failed_kmer_max_entry_count - _failed_kmer_entries_count > no_of_keys / 2) {
			// there is a possibility that all keys would not get inserted.
			return false;
		}

		uint32_t no_of_keys_per_stream = (no_of_keys + _no_of_streams - 1) / _no_of_streams;
		for (uint32_t i = 0; i < _no_of_streams; i++) {
			InsertToHashTable(d_input, no_of_keys_per_stream, _cuda_streams[i], _d_kmer_db, _kmer_db_max_record_count, _d_cuda_counters + (i * 2));
		}

		UpdateCounters();

		return true;
	}

	void Extract(std::vector<KmerKeyValue<key_size> >& kmer_entries) const {
		// set device id
		CUDA_CHECK_RETURN(cudaSetDevice(_device_id));

		uint64_t read_size = 0;
		uint64_t temp_buffer_size = sizeof(KmerKeyValue<key_size>) * 1000;
		char* temp = new char[temp_buffer_size];
		for (uint64_t index = 0; index < _kmer_db_size; index += read_size) {
			if (index + temp_buffer_size < _kmer_db_size) {
				CUDA_CHECK_RETURN(cudaMemcpy(temp, _d_kmer_db + index, temp_buffer_size, cudaMemcpyDeviceToHost));
				read_size = temp_buffer_size;
			} else {
				CUDA_CHECK_RETURN(cudaMemcpy(temp, _d_kmer_db + index, _kmer_db_size - index, cudaMemcpyDeviceToHost));
			    read_size = _kmer_db_size - index;
			}

			for (int64_t j = 0; j < read_size; j += sizeof(KmerKeyValue<key_size>)) {
				if (((KmerKeyValue<key_size>*)(temp + j))->getCount() > 0) {
					kmer_entries.push_back(*(KmerKeyValue<key_size>*)(temp + j));
				}
			}
		}

		// TODO : Implement kmer extraction logic
		// need to copy the kmer db and failed keys both and clear the table.
	}

	void UpdateCounters() {
		CUDA_CHECK_RETURN(cudaDeviceSynchronize());
		_kmer_db_record_count = 0;
		_failed_kmer_entries_count = 0;

		CUDA_CHECK_RETURN(cudaMemcpy(_h_cuda_counters, _d_cuda_counters, _cuda_counters_per_stream * _no_of_streams * sizeof(uint64_t), cudaMemcpyDeviceToHost));
		for (uint32_t i = 0; i < _cuda_counters_per_stream * _no_of_streams; i += _cuda_counters_per_stream) {
			_kmer_db_record_count += _h_cuda_counters[i + 0];
			_failed_kmer_entries_count += _h_cuda_counters[i + 1];
		}

		for (uint32_t i = 0; i < _cuda_counters_per_stream * _no_of_streams; i++) {
			printf("========================%"PRIu64"\n", _h_cuda_counters[i]);
		}
		printf("===========s:%"PRIu64", f:%"PRIu64"\n", _kmer_db_record_count, _failed_kmer_entries_count);
	}

	uint64_t GetSuccessCount() {
		return _kmer_db_record_count;
	}

	uint64_t GetFailureCount() {
		return _failed_kmer_entries_count;
	}

	void Dump() {

		std::ofstream output_file_0("/tmp/1/out.dump");


		// set device id
		CUDA_CHECK_RETURN(cudaSetDevice(_device_id));

		CUDA_CHECK_RETURN(cudaDeviceSynchronize());

		uint64_t read_size = 0;
		uint64_t temp_buffer_size = sizeof(KmerKeyValue<key_size>) * 1000;
		char* temp = new char[temp_buffer_size];
		for (uint64_t index = 0; index < _kmer_db_size; index += read_size) {
			if (index + temp_buffer_size < _kmer_db_size) {
				(cudaMemcpy(temp, ((char*)_d_kmer_db) + index, temp_buffer_size, cudaMemcpyDeviceToHost));
				read_size = temp_buffer_size;
			} else {
				CUDA_CHECK_RETURN(cudaMemcpy(temp, ((char*)_d_kmer_db) + index, _kmer_db_size - index, cudaMemcpyDeviceToHost));
				read_size = _kmer_db_size - index;
			}

			for (int64_t j = 0; j < read_size; j += sizeof(KmerKeyValue<key_size>)) {
				if (((KmerKeyValue<key_size>*)(temp + j))->getCount() > 0) {
					output_file_0.write((temp + j), sizeof(KmerKeyValue<key_size>) - sizeof(KmerKey<key_size>));
					//output_file_0.write(&((uint32_t)((KmerKeyValue<key_size>*)(temp + j))->getCount()), 4);
					uint32_t count = (uint32_t)((KmerKeyValue<key_size>*)(temp + j))->getCount();
					output_file_0.write((char*)&count, sizeof(uint32_t));
				}
			}
		}

		output_file_0.close();

	}

	void TempDump() {
		std::ofstream output_file_0("/tmp/1/0", std::ofstream::out | std::ofstream::app);
		std::ofstream output_file_1("/tmp/1/1", std::ofstream::out | std::ofstream::app);
		std::ofstream output_file_2("/tmp/1/2", std::ofstream::out | std::ofstream::app);
		std::ofstream output_file_3("/tmp/1/3", std::ofstream::out | std::ofstream::app);

		char* data = new char[_kmer_db_size];
		CUDA_CHECK_RETURN(cudaMemcpy(data, ((char*)_d_kmer_db), _kmer_db_size, cudaMemcpyDeviceToHost));

		std::ofstream* output_file = NULL;
		for (int64_t j = 0; j < _kmer_db_size; j += sizeof(KmerKeyValue<key_size>)) {
			KmerKeyValue<key_size>* key_value = (KmerKeyValue<key_size>*)(data + j);
			if (key_value->getCount() > 0) {
				int32_t first_letter = *(uint64_t*)(&key_value->getKey()) >> 62;
				switch(first_letter) {
				case 0:
					output_file = &output_file_0;
					break;
				case 1:
					output_file = &output_file_1;
					break;
				case 2:
					output_file = &output_file_2;
					break;
				case 3:
					output_file = &output_file_3;
					break;
				default:
					printf("ERROR: Wring starting letter : %i", first_letter);
					continue;
				}
				output_file->write((data + j), sizeof(KmerKeyValue<key_size>) - sizeof(KmerKey<key_size>));
				uint32_t count = (uint32_t)key_value->getCount();
				output_file->write((char*)&count, sizeof(uint32_t));
			}
		}

		delete[] data;

		output_file_0.close();
		output_file_1.close();
		output_file_2.close();
		output_file_3.close();
	}
};

#endif /* COUNTINGHASHTABLEBASE_H_ */
