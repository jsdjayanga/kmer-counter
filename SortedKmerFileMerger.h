/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   SortedKmerFileMerger.h
 * Author: jayanga
 *
 * Created on December 8, 2016, 6:26 PM
 */

#ifndef SORTEDKMERFILEMERGER_H
#define SORTEDKMERFILEMERGER_H

#include <string>
#include <list>

#include "SortedKmerFile.h"

class SortedKmerFileMerger {
public:
    SortedKmerFileMerger(std::string output_filename);
    SortedKmerFileMerger(const SortedKmerFileMerger& orig);
    virtual ~SortedKmerFileMerger();
    void Merge(std::list<std::string> files, uint32_t kmer_length);
private:
    SortedKmerFile* GetLowest(std::list<SortedKmerFile*>& sorted_kmer_files);
    std::string _output_filename;
};

#endif /* SORTEDKMERFILEMERGER_H */

