/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   GPUHandler.h
 * Author: jayanga
 *
 * Created on October 21, 2016, 10:22 PM
 */

#pragma once

#include <stdint.h>

using namespace std;

int64_t processKMers(const char* input, int64_t kmerLength, int64_t inputSize,
		int64_t lineLength);

