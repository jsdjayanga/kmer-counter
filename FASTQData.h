/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FASTQData.h
 * Author: jayanga
 *
 * Created on October 19, 2016, 5:48 PM
 */

#pragma once

#include <stdint.h>
#include <fstream>

using namespace std;

class FASTQData {
public:
    FASTQData(char* data, int64_t size, int64_t lineLength);
    FASTQData(const FASTQData& orig);
    virtual ~FASTQData();
    
    //void read();
    const char* getData();
    void setLineLength(int64_t _lineLength);
    int64_t getLineLength() const;
    void setSize(int64_t _size);
    int64_t getSize() const;
    void setData(char* _data);
private:
    char*   _data;
    int64_t _size;
    int64_t _lineLength;
};


