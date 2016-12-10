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

#include "SortedKmerMerger.h"

using namespace std;

SortedKmerMerger::SortedKmerMerger(uint32_t kmer_length) {
	_kmer_length = kmer_length;

	_kmer_store_size = _kmer_length / 32;
	if (_kmer_length % 32 > 0) {
		_kmer_store_size++;
	}
	_kmer_store_size *= 8;
	_kmer_store_size += 8;
}

SortedKmerMerger::SortedKmerMerger(const SortedKmerMerger& orig) {
}

SortedKmerMerger::~SortedKmerMerger() {
}

void SortedKmerMerger::Merge(std::list<std::pair<char*, uint64_t> > kmer_list, std::string filename) {
    list<SortedKmerArray*> sorted_kmer_arrays;

    for (std::list<std::pair<char*, uint64_t> >::iterator it = kmer_list.begin(); it != kmer_list.end(); ++it) {
        sorted_kmer_arrays.push_back(new SortedKmerArray(it->first, it->second));
    }

    ofstream output_file(filename.c_str());
    char* data = GetLowest(sorted_kmer_arrays);
    while (data != NULL) {
        output_file.write(data, _kmer_store_size - sizeof(uint64_t));
        uint32_t count = (uint32_t)*(uint64_t*)(data + _kmer_store_size - sizeof(uint64_t));
        output_file.write((char*)&count, sizeof(uint32_t));
        
//        cout << "F" << ":" << *(uint64_t*)data << "|" << count << endl;
        
        data = GetLowest(sorted_kmer_arrays);
    }
    output_file.close();

    for (std::list<std::pair<char*, uint64_t> >::iterator it = kmer_list.begin(); it != kmer_list.end(); ++it) {
		delete[] it->first;
	}
}

char* SortedKmerMerger::GetLowest(std::list<SortedKmerArray*>& sorted_kmer_arrays) {
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
        
        if (lessThan(current->_data + current->_index, lowest->_data + lowest->_index, _kmer_length)) {
            lowest = current;
            lowest_ite = ite;
        } else if (equals(current->_data + current->_index, lowest->_data + lowest->_index, _kmer_length)) {
            *(uint64_t*) (lowest->_data + lowest->_index + _kmer_store_size - sizeof(uint64_t)) =
            		*(uint64_t*) (lowest->_data + lowest->_index + _kmer_store_size - sizeof(uint64_t)) +
                    *(uint64_t*) (current->_data + current->_index + _kmer_store_size - sizeof(uint64_t));
            current->_index += _kmer_store_size;
            if (current->_index >= current->_length) {
                list<SortedKmerArray*>::iterator ite_to_delete = ite;
                ite++;
                sorted_kmer_arrays.erase(ite_to_delete);
                continue;
            }
        }
        ite++;
    }

    uint64_t lowest_index = lowest->_index;
    lowest->_index += _kmer_store_size;
    if (lowest->_index >= lowest->_length) {
        sorted_kmer_arrays.erase(lowest_ite);
    }
    return (char*) lowest->_data + lowest_index;
}
