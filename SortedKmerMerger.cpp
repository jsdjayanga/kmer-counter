/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   SortedKmerMerger.cpp
 * Author: jayanga
 * 
 * Created on December 5, 2016, 9:26 PM
 */

#include <fstream>
#include <iostream>
#include "KMerCounterUtils.h"
#include <cstring>

#include "SortedKmerMerger.h"

using namespace std;

SortedKmerMerger::SortedKmerMerger(uint32_t kmer_length) {
	_kmer_length = kmer_length;

	_key_size_in_longs = _kmer_length / 32;
	if (_kmer_length % 32 > 0) {
		_key_size_in_longs++;
	}
	_kmer_store_size = _key_size_in_longs * 8;
	_kmer_store_size += 8;

	_buffer_index = 0;
	_buffer_size = (uint64_t)1 * 1024 * 1024 * 1024;
	_buffer = new char[_buffer_size];
}

SortedKmerMerger::SortedKmerMerger(const SortedKmerMerger& orig) {
}

SortedKmerMerger::~SortedKmerMerger() {
	delete[] _buffer;
}

void SortedKmerMerger::Merge(std::list<std::pair<char*, uint64_t> > kmer_list, std::string filename) {
    list<SortedKmerArray*> sorted_kmer_arrays;
    list<SortedKmerArray*> ararys_to_delete;

    for (std::list<std::pair<char*, uint64_t> >::iterator it = kmer_list.begin(); it != kmer_list.end(); ++it) {
        sorted_kmer_arrays.push_back(new SortedKmerArray(it->first, it->second));
    }

    ofstream output_file(filename.c_str());
    char* data = GetLowest(sorted_kmer_arrays, ararys_to_delete);
    uint64_t result_kmer_store_size = _kmer_store_size - sizeof(uint64_t) + sizeof(uint32_t);
    while (data != NULL) {

    	if (_buffer_index + _kmer_store_size > _buffer_size) {
			output_file.write(_buffer, _buffer_index);
			_buffer_index = 0;
		}

    	memcpy(_buffer + _buffer_index, data, _kmer_store_size - sizeof(uint64_t));
		uint32_t count = (uint32_t)*(uint64_t*)(data + _kmer_store_size - sizeof(uint64_t));
		memcpy(_buffer + _buffer_index + _kmer_store_size - sizeof(uint64_t), (char*)&count, sizeof(uint32_t));
		_buffer_index += result_kmer_store_size;

//        output_file.write(data, _kmer_store_size - sizeof(uint64_t));
//        uint32_t count = (uint32_t)*(uint64_t*)(data + _kmer_store_size - sizeof(uint64_t));
//        output_file.write((char*)&count, sizeof(uint32_t));
        
//        cout << "F" << ":" << *(uint64_t*)data << "|" << count << endl;
        
        data = GetLowest(sorted_kmer_arrays, ararys_to_delete);
    }

    output_file.write(_buffer, _buffer_index);
    _buffer_index = 0;

    output_file.close();

    for (list<SortedKmerArray*>::iterator it = ararys_to_delete.begin(); it != ararys_to_delete.end(); ++it) {
		delete *it;
	}
}

char* SortedKmerMerger::GetLowest(std::list<SortedKmerArray*>& sorted_kmer_arrays, list<SortedKmerArray*>& ararys_to_delete) {
    if (sorted_kmer_arrays.empty()) {
        return NULL;
    }

    list<SortedKmerArray*>::iterator ite = sorted_kmer_arrays.begin();
    list<SortedKmerArray*>::iterator lowest_ite = ite;
    SortedKmerArray* lowest = *ite;
    
    for (ite++; ite != sorted_kmer_arrays.end();) {
        SortedKmerArray* current = *ite;
        
//        cout << *(uint64_t*)(lowest->_data + lowest->_index) << "|" << *(uint64_t*)(lowest->_data + lowest->_index + 8)
//                << "|" << *(uint64_t*)(current->_data + current->_index) << endl;
        
        if (lessThan(current->_data + current->_index, lowest->_data + lowest->_index, _key_size_in_longs)) {
            lowest = current;
            lowest_ite = ite;
        } else if (equals(current->_data + current->_index, lowest->_data + lowest->_index, _key_size_in_longs)) {
            *(uint64_t*) (lowest->_data + lowest->_index + _kmer_store_size - sizeof(uint64_t)) =
            		*(uint64_t*) (lowest->_data + lowest->_index + _kmer_store_size - sizeof(uint64_t)) +
                    *(uint64_t*) (current->_data + current->_index + _kmer_store_size - sizeof(uint64_t));
            current->_index += _kmer_store_size;
            if (current->_index >= current->_length) {
                list<SortedKmerArray*>::iterator ite_to_delete = ite;
                ite++;
                ararys_to_delete.push_back(*ite_to_delete);
                sorted_kmer_arrays.erase(ite_to_delete);
                continue;
            }
        }
        ite++;
    }

    uint64_t lowest_index = lowest->_index;
    lowest->_index += _kmer_store_size;
    char* data = (char*) lowest->_data + lowest_index;
    if (lowest->_index >= lowest->_length) {
    	ararys_to_delete.push_back(*lowest_ite);
        sorted_kmer_arrays.erase(lowest_ite);
    }
    return data;
}
