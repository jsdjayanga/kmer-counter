/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FASTQFileReader.cpp
 * Author: jayanga
 * 
 * Created on October 19, 2016, 5:45 PM
 */

#include <iostream>
#include <cstring>
#include "FASTQFileReader.h"

FASTQFileReader::FASTQFileReader(string filename, int64_t fileSize) {
    this->_filename = filename;
    this->_fileSize = fileSize;

    ifstream fileStreamB;
    fileStreamB.open(_filename.c_str(), ios::ate | ios::binary);
    if (fileStreamB.is_open()) {
        _fileSize = fileStreamB.tellg();
    }

    _fileStream.open(_filename.c_str());
    if (_fileStream.is_open()) {
        string line1;
        string line2;
        getline(_fileStream, line1);
        getline(_fileStream, line2);

        _lineLength = line2.size();

        _fileStream.seekg(0);
    }

    cout << _fileSize << "|" << _lineLength << endl;
}

FASTQFileReader::FASTQFileReader(const FASTQFileReader& orig) {
}

FASTQFileReader::~FASTQFileReader() {
}

FASTQData* FASTQFileReader::readData(int64_t chunkSize) {
    char* data = new char[chunkSize];

    int64_t chunk_offset = 0;
    string temp;
    string line;
    getline(_fileStream, temp);
    getline(_fileStream, line);
    while (line.length() != 0 && chunk_offset + temp.length() < chunkSize) {
        const char* cline = line.c_str();
        if (cline != NULL && line.length() > 0 && cline[0] == '+') {
            const char* cline = temp.c_str();

            if (chunk_offset + temp.length() < chunkSize) {
                strncpy(&data[chunk_offset], cline, temp.length());
                chunk_offset += temp.length();

//                if (chunk_offset < 200) {
//                    cout << temp << endl;
//                }

                getline(_fileStream, temp);
                getline(_fileStream, line);
            } else {
                _fileStream.seekg(_fileStream.tellg() - temp.length());
            }
        } else {
            temp = line;
            getline(_fileStream, line);
        }

        //cout << temp << "======" << line << "====" << chunk_offset << endl;

    }
    data[chunk_offset] = '\0';
    // cout << data << endl;
    //cout << temp << endl;

    FASTQData* fastqData = new FASTQData(data, chunk_offset, _lineLength);
    return fastqData;
}

bool FASTQFileReader::isComplete() {
    return _fileStream.tellg() + _lineLength > _fileSize;
}

int64_t FASTQFileReader::getLineLength() {
    return _lineLength;
}
