/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   SortedKmerMerger.h
 * Author: jayanga
 *
 * Created on December 5, 2016, 9:26 PM
 */

#ifndef SORTEDKMERMERGER_H
#define SORTEDKMERMERGER_H

#include <string>
#include <stdint.h>
#include <list>

struct SortedKmerArray {

    SortedKmerArray(char* data, uint64_t length) {
        _data = data;
        _index = 0;
        _length = length;
    }
    char* _data;
    uint64_t _index;
    uint64_t _length;
};

class SortedKmerMerger {
public:
    SortedKmerMerger();
    SortedKmerMerger(const SortedKmerMerger& orig);
    virtual ~SortedKmerMerger();
    void Merge(std::list<std::pair<char*, uint64_t> > kmer_list, std::string filename);
private:
    char* GetLowest(std::list<SortedKmerArray*>& sorted_kmer_arrays);
};

#endif /* SORTEDKMERMERGER_H */

