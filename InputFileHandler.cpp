/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   InputFileHandler.cpp
 * Author: jayanga
 * 
 * Created on October 19, 2016, 8:20 PM
 */

#include "InputFileHandler.h"
#include <dirent.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string.h>

InputFileHandler::InputFileHandler(string fileDirectory) {
    this->_fileDirectory = fileDirectory;

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

                FASTQFileReader* fastqFileReader = new FASTQFileReader(filename, temp.tellg());
                _fileReaders.push_back(fastqFileReader);
            }
        }
        closedir(dir);
    } else {
        cout << "Couldn't open directory : " << _fileDirectory << endl;
    }
}

InputFileHandler::InputFileHandler(const InputFileHandler& orig) {
}

InputFileHandler::~InputFileHandler() {
}

//list<InputFileDetails*>& InputFileHandler::getFileList() {
//    DIR *dir;
//    struct dirent *ent;
//    if ((dir = opendir(_fileDirectory.c_str())) != NULL) {
//        /* print all the files and directories within directory */
//        while ((ent = readdir(dir)) != NULL) {
//            if (strncmp(".", ent->d_name, 1) != 0 && strncmp("..", ent->d_name, 2) != 0) {
//                //printf("%s\n", ent->d_name);
//                
//                ifstream temp;
//                string filename = _fileDirectory + "/" + ent->d_name;
//                temp.open(filename.c_str(), ios::ate | ios::binary);
//
//                InputFileDetails* fileDetails = new InputFileDetails(ent->d_name, temp.tellg());
//                _inputFileDetails.push_back(fileDetails);
//            }
//        }
//        closedir(dir);
//    } else {
//        cout << "Couldn't open directory : " << _fileDirectory << endl;
//    }
//    return _inputFileDetails;
//}

FASTQData* InputFileHandler::read(int64_t chunkSize) {
    if (!_fileReaders.empty()) {
        FASTQFileReader* fastqFileReader = _fileReaders.front();
        if (fastqFileReader != NULL) {
            FASTQData* fastqData = fastqFileReader->readData(chunkSize);
            if (fastqFileReader->isComplete() || fastqData->getSize() == 0) {
                _fileReaders.pop_front();
            }
            return fastqData;
        }
    } else {
        return NULL;
    }
}

int64_t InputFileHandler::getLineLength() {
    if (!_fileReaders.empty()) {
        FASTQFileReader* fastqFileReader = _fileReaders.front();
        if (fastqFileReader != NULL) {
            return fastqFileReader->getLineLength();
        }
    }
    return 0;
}
