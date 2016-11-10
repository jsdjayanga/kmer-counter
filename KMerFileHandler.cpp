/*
 * KMerFileHandler.cpp
 *
 *  Created on: Nov 10, 2016
 *      Author: jayanga
 */

#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include "KMerFileHandler.h"

KMerFileHandler::KMerFileHandler(string fileDirectory, uint64_t recordLength, uint64_t recordCount) {
	// TODO Auto-generated constructor stub
	this->_fileDirectory = fileDirectory;
	this->_recordLength = recordLength;
	this->_recordCount = recordCount;
}

KMerFileHandler::~KMerFileHandler() {
	// TODO Auto-generated destructor stub
}

list<KMerFileReader*>& KMerFileHandler::getKmerFileList() {
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir(_fileDirectory.c_str())) != NULL) {
		/* print all the files and directories within directory */
		while ((ent = readdir(dir)) != NULL) {
			if (strncmp(".", ent->d_name, 1) != 0 && strncmp("..", ent->d_name, 2) != 0) {
				//printf("%s\n", ent->d_name);

				ifstream temp;
				string filename = _fileDirectory + "/" + ent->d_name;
				temp.open(filename.c_str(), ios::ate | ios::binary);

				KMerFileReader* kMerFileReader = new KMerFileReader(filename, temp.tellg(), _recordLength, _recordCount);
				_kmerFileReaders.push_back(kMerFileReader);
			}
		}
		closedir(dir);
	} else {
		cout << "Couldn't open directory : " << _fileDirectory << endl;
	}

	return _kmerFileReaders;
}
