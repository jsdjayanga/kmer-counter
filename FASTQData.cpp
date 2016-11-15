/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FASTQData.cpp
 * Author: jayanga
 * 
 * Created on October 19, 2016, 5:48 PM
 */

#include "FASTQData.h"

FASTQData::FASTQData(char* data, int64_t size, int64_t lineLength) {
    _data = data;
    _size = size;
    _lineLength = lineLength;
}

FASTQData::FASTQData(const FASTQData& orig) {
}

FASTQData::~FASTQData() {
	delete[] _data;
}

//void FASTQData::read() {
//    int64_t offset = 0;
//    string temp;
//    string line2;
////    getline(_fileStream, temp);
////    getline(_fileStream, line2);
////    while (offset < _size && line2.length() != 0) {
////        const char* cline2 = line2.c_str();
////        if (cline2 != NULL && line2.length() > 0 && cline2[0] == '+') {
////            const char* cline = temp.c_str();
////
////            strncpy(&_data[offset], cline, temp.length());
////            offset += temp.length();
////            //cout << temp << endl;
////
////            getline(_fileStream, temp);
////            getline(_fileStream, line2);
////        } else {
////            temp = line2;
////            getline(_fileStream, line2);
////        }
////
////        //cout << temp << "======" << line2 << "====" << offset << endl;
////
////    }
////    _data[offset] = '\0';
////    cout << _data << endl;
//}

const char* FASTQData::getData() {
    return _data;
}

void FASTQData::setLineLength(int64_t _lineLength) {
    this->_lineLength = _lineLength;
}

int64_t FASTQData::getLineLength() const {
    return _lineLength;
}

void FASTQData::setSize(int64_t _size) {
    this->_size = _size;
}

int64_t FASTQData::getSize() const {
    return _size;
}

void FASTQData::setData(char* _data) {
    this->_data = _data;
}

