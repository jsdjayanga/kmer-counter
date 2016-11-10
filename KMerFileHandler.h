/*
 * KMerFileHandler.h
 *
 *  Created on: Nov 10, 2016
 *      Author: jayanga
 */

#pragma once

#include <list>
#include "KMerFileReader.h"

using namespace std;

class KMerFileHandler {
public:
	KMerFileHandler(string fileDirectory, uint64_t recordLength, uint64_t recordCount);
	virtual ~KMerFileHandler();
	list<KMerFileReader*>& getKmerFileList();
private:
	string _fileDirectory;
	list<KMerFileReader*> _kmerFileReaders;
};

