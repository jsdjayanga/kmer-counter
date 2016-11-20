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
#include "KMerPrinter.h"

using namespace std;

Options* getOptions(int argc, char** argv) {
	Options* options = new Options();
	options->SetInputFileDirectory("/home/jayangad/data/1");
	options->SetGpuMemoryLimit(100000000);
	options->setTempFileLocation("/tmp/1");
	options->setOutputFile("/tmp/2/output.bin");

	for (int i = 0; i < argc; i++) {
		if (strncmp(argv[i], "kmerLength=", 11) == 0) {
			options->SetKmerLength(strtoull(argv[i] + 11, NULL, 10));
			cout << "Updating KmerLength=" << options->GetKmerLength() << endl;
		}

		if (strncmp(argv[i], "gpuMemoryLimit=", 15) == 0) {
			options->SetGpuMemoryLimit(strtoull(argv[i] + 15, NULL, 10));
			cout << "Updating Gpu Memory Limit=" << options->GetGpuMemoryLimit() << endl;
		}

		if (strncmp(argv[i], "inputFileLocation=", 18) == 0) {
			options->SetInputFileDirectory(argv[i] + 18);
			cout << "Updating Input File Location='" << options->GetInputFileDirectory() << "'" << endl;
		}

		if (strncmp(argv[i], "tempFileLocation=", 17) == 0) {
			options->setTempFileLocation(argv[i] + 17);
			cout << "Updating Temp File Location='" << options->getTempFileLocation() << "'" << endl;
		}

		if (strncmp(argv[i], "outputFile=", 11) == 0) {
			options->setOutputFile(argv[i] + 11);
			cout << "Updating Output File='" << options->getOutputFile() << "'" << endl;
		}

		if (strncmp(argv[i], "noOfMergersAtOnce=", 18) == 0) {
			options->setNoOfMergersAtOnce(atoi(argv[i] + 18));
			cout << "Updating No Of Mergers At Once='" << options->getNoOfMergersAtOnce() << "'" << endl;
		}

		if (strncmp(argv[i], "noOfMergeThreads=", 17) == 0) {
					options->setNoOfMergeThreads(atoi(argv[i] + 17));
					cout << "Updating No Of Merge Threads='" << options->getNoOfMergeThreads() << "'" << endl;
				}
	}

	return options;
}

/*
 * 
 */
int main(int argc, char** argv) {
	cout << "### kmer-counter application ###" << endl;

	if (argc == 5 && strncmp(argv[1], "print", 5) == 0) {
		KMerPrinter* kMerPrinter = new KMerPrinter(argv[2], argv[3], atoll(argv[4]));
		kMerPrinter->print();
		return 0;
	}

	Options* options = getOptions(argc, argv);
	KMerCounter* kmerCounter = new KMerCounter(options);
	kmerCounter->Start();

	return 0;
}

