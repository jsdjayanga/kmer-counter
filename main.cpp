/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: jayanga
 *
 * Created on October 19, 2016, 5:40 PM
 */

#include <cstdlib>
#include <iostream>
#include <string.h>

#include "FASTQFileReader.h"
#include "Options.h"
#include "KMerCounter.h"

using namespace std;

Options* getOptions(int argc, char** argv) {
	Options* options = new Options();
	options->SetInputFileDirectory("/home/jayangad/data/1");
	options->SetGpuMemoryLimit(100000000);

	for (int i = 0; i < argc; i++) {
		if (strncmp(argv[i], "kmerLength=", 11) == 0) {
			options->SetKmerLength(strtoull(argv[i] + 11, NULL, 10));
			cout << "Updating KmerLength=" << options->GetKmerLength() << endl;
		}

		if (strncmp(argv[i], "gpuMemoryLimit=", 15) == 0) {
			options->SetGpuMemoryLimit(strtoull(argv[i] + 15, NULL, 10));
			cout << "Updating Gpu Memory Limit=" << options->GetGpuMemoryLimit() << endl;
		}
	}

	return options;
}

/*
 * 
 */
int main(int argc, char** argv) {
	cout << "### kmer-counter application ###" << endl;

	Options* options = getOptions(argc, argv);
	KMerCounter* kmerCounter = new KMerCounter(options);
	kmerCounter->Start();

	return 0;
}

