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
#include "SortedKmerFileMerger.h"
#include "SortedKmerFile.h"

using namespace std;

SortedKmerFileMerger::SortedKmerFileMerger(string output_filename) {
    _output_filename = output_filename;
}

SortedKmerFileMerger::SortedKmerFileMerger(const SortedKmerFileMerger& orig) {
}

SortedKmerFileMerger::~SortedKmerFileMerger() {
}

void SortedKmerFileMerger::Merge(list<string> files, uint32_t kmer_length) {
    list<SortedKmerFile*> sorted_kmer_files;
    
    for (list<string>::iterator ite = files.begin(); ite != files.end(); ite++) {
        SortedKmerFile* sf = new SortedKmerFile(*ite, kmer_length);
        if (sf->Peek() != NULL) {
            sorted_kmer_files.push_back(sf);
        }
    }
    
    ofstream output_file(_output_filename.c_str());
    SortedKmerFile* lowest = GetLowest(sorted_kmer_files);
    while (lowest != NULL) {
        output_file.write(lowest->Peek(), sizeof(uint64_t));
        uint32_t count = (uint32_t)*(uint64_t*)(lowest->Peek() + sizeof(uint64_t));
        output_file.write((char*)&count, sizeof(uint32_t));        
        
//        cout << *(uint64_t*)lowest->Peek() << " " << count << endl;
        
        lowest->Pop();
         
        lowest = GetLowest(sorted_kmer_files);
    }
    output_file.close();
}

SortedKmerFile* SortedKmerFileMerger::GetLowest(std::list<SortedKmerFile*>& sorted_kmer_files) {
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
            sorted_kmer_files.erase(ite_to_delete);
            lowest = current;
            continue;
        }
        
        if (current->Peek() == NULL) {
            list<SortedKmerFile*>::iterator ite_to_delete = ite;
            ite++;
            sorted_kmer_files.erase(ite_to_delete);
            continue;
        }
        
//        cout << *(uint64_t*)(lowest->_data + lowest->_index) << "|" << *(uint64_t*)(lowest->_data + lowest->_index + 8)
//                << "|" << *(uint64_t*)(current->_data + current->_index) << endl;
        
        if (*(uint64_t*) (current->Peek()) < *(uint64_t*) (lowest->Peek())) {
            lowest = current;
            lowest_ite = ite;
        } else if (*(uint64_t*) (current->Peek()) == *(uint64_t*) (lowest->Peek())) {
            lowest->UpdatePeek(*(uint32_t*)(current->Peek() + sizeof(uint64_t)));
            current->Pop();
//            *(uint64_t*) (lowest->_data + lowest->_index + sizeof (uint64_t)) = *(uint64_t*) (lowest->_data + lowest->_index + sizeof (uint64_t)) +
//                    *(uint64_t*) (current->_data + current->_index + sizeof (uint64_t));
//            current->_index += 2 * sizeof (uint64_t);
            if (current->Peek() == NULL) {
                list<SortedKmerFile*>::iterator ite_to_delete = ite;
                ite++;
                sorted_kmer_files.erase(ite_to_delete);
                continue;
            }
        }
        ite++;
    }

    if (lowest->Peek() != NULL) {
        return lowest;
    }
    return NULL;
}
