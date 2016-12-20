/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   SortedKmerFileMerger.cpp
 * Author: jayanga
 * 
 * Created on December 8, 2016, 6:26 PM
 */

#include <iostream>
#include <cstring>
#include "SortedKmerFileMerger.h"
#include "SortedKmerFile.h"
#include "KMerCounterUtils.h"

using namespace std;

SortedKmerFileMerger::SortedKmerFileMerger(string output_filename, uint32_t kmer_length) {
    _output_filename = output_filename;
    
    _buffer_index = 0;
    _buffer_size = (uint64_t)1 * 1024 * 1024 * 1024;
    _buffer = new char[_buffer_size];

    _key_size_in_longs = kmer_length / 32;
	if (kmer_length % 32 > 0) {
		_key_size_in_longs++;
	}
}

SortedKmerFileMerger::SortedKmerFileMerger(const SortedKmerFileMerger& orig) {
}

SortedKmerFileMerger::~SortedKmerFileMerger() {
	delete[] _buffer;
}

void SortedKmerFileMerger::Merge(list<string> files, uint32_t kmer_length) {
    list<SortedKmerFile*> sorted_kmer_files;
    list<SortedKmerFile*> files_to_delete;
    
    for (list<string>::iterator ite = files.begin(); ite != files.end(); ite++) {
        SortedKmerFile* sf = new SortedKmerFile(*ite, kmer_length);
        if (sf->Peek() != NULL) {
            sorted_kmer_files.push_back(sf);
        } else {
        	delete sf;
        }
    }
    
    uint32_t kmer_store_size = kmer_length / 32;
	if (kmer_length % 32 > 0) {
		kmer_store_size++;
	}
	kmer_store_size *= 8;
	kmer_store_size += 4; // This should be 4: it is read from the file

    ofstream output_file(_output_filename.c_str());
    SortedKmerFile* lowest = GetLowest(sorted_kmer_files, files_to_delete, kmer_length, kmer_store_size);
    while (lowest != NULL) {
        
        if (_buffer_index + kmer_store_size > _buffer_size) {
            output_file.write(_buffer, _buffer_index);
            _buffer_index = 0;
        }
        
        memcpy(_buffer + _buffer_index, lowest->Peek(), kmer_store_size - sizeof(uint32_t));
        uint32_t count = (uint32_t)*(uint64_t*)(lowest->Peek() + kmer_store_size - sizeof(uint32_t));
        memcpy(_buffer + _buffer_index + kmer_store_size - sizeof(uint32_t), (char*)&count, sizeof(uint32_t));
        _buffer_index += kmer_store_size;
        
//        output_file.write(lowest->Peek(), sizeof(uint64_t));
//        uint32_t count = (uint32_t)*(uint64_t*)(lowest->Peek() + sizeof(uint64_t));
//        output_file.write((char*)&count, sizeof(uint32_t));        
        
//        cout << *(uint64_t*)lowest->Peek() << " " << count << endl;
        
        lowest->Pop();
         
        lowest = GetLowest(sorted_kmer_files, files_to_delete, kmer_length, kmer_store_size);
    }
    
	output_file.write(_buffer, _buffer_index);
	_buffer_index = 0;

    output_file.close();

    for (list<SortedKmerFile*>::iterator it = files_to_delete.begin(); it != files_to_delete.end(); ++it) {
		delete *it;
	}
}

SortedKmerFile* SortedKmerFileMerger::GetLowest(std::list<SortedKmerFile*>& sorted_kmer_files, list<SortedKmerFile*>& files_to_delete, uint32_t kmer_length, uint32_t kmer_store_size) {
    if (sorted_kmer_files.empty()) {
        return NULL;
    }

    list<SortedKmerFile*>::iterator ite = sorted_kmer_files.begin();
    list<SortedKmerFile*>::iterator lowest_ite = ite;
    SortedKmerFile* lowest = *ite;
    
    for (ite++; ite != sorted_kmer_files.end();) {
        SortedKmerFile* current = *ite;
        
        if (lowest->Peek() == NULL) {
            list<SortedKmerFile*>::iterator ite_to_delete = ite;
            ite++;
            files_to_delete.push_back(*ite_to_delete);
            sorted_kmer_files.erase(ite_to_delete);
            lowest = current;
            continue;
        }
        
        if (current->Peek() == NULL) {
            list<SortedKmerFile*>::iterator ite_to_delete = ite;
            ite++;
            files_to_delete.push_back(*ite_to_delete);
            sorted_kmer_files.erase(ite_to_delete);
            continue;
        }
        
//        cout << *(uint64_t*)(lowest->_data + lowest->_index) << "|" << *(uint64_t*)(lowest->_data + lowest->_index + 8)
//                << "|" << *(uint64_t*)(current->_data + current->_index) << endl;
        
        if (lessThan(current->Peek(), lowest->Peek(), _key_size_in_longs)) {
            lowest = current;
            lowest_ite = ite;
        } else if (equals(current->Peek(), lowest->Peek(), _key_size_in_longs)) {
            lowest->UpdatePeek(*(uint32_t*)(current->Peek() + kmer_store_size - sizeof(uint32_t)));
            current->Pop();
//            *(uint64_t*) (lowest->_data + lowest->_index + sizeof (uint64_t)) = *(uint64_t*) (lowest->_data + lowest->_index + sizeof (uint64_t)) +
//                    *(uint64_t*) (current->_data + current->_index + sizeof (uint64_t));
//            current->_index += 2 * sizeof (uint64_t);
            if (current->Peek() == NULL) {
                list<SortedKmerFile*>::iterator ite_to_delete = ite;
                ite++;
                files_to_delete.push_back(*ite_to_delete);
                sorted_kmer_files.erase(ite_to_delete);
                continue;
            }
        }
        ite++;
    }

    if (lowest->Peek() != NULL) {
        return lowest;
    }

    files_to_delete.push_back(lowest);
    return NULL;
}
