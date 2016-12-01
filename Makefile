## -*- Makefile -*-
##
## User: jayanga
## Time: Oct 21, 2016 10:08:48 PM
## Makefile created by Oracle Developer Studio.
##
## This file is generated automatically.
##


all: kmer-counter

## Target: kmer-counter
kmer-counter: cuda
	g++ -ggdb -O0 -std=c++11 -Wl,--no-as-needed -I /usr/local/cuda/include/ -I /home/jayangad/source/sparsehash/build/include/ -I /home/jayangad/source/tbb2017_20161004oss/include -o kmer-counter main.cpp FASTQData.cpp FASTQFileReader.cpp FileDump.cpp InputFileHandler.cpp KMerCounter.cpp Options.cpp KMerFileHandler.cpp KMerFileReader.cpp KMerFileMerger.cpp KMerPrinter.cpp SortedKMerFile.cpp KMerFileMergeHandler.cpp GPUHandler.o CountingHashTable.o -L/usr/local/cuda/lib64 -L/home/jayangad/source/tbb2017_20161004oss/lib/intel64/gcc4.7/ -lcudart -pthread -ltbb

cuda:
	nvcc -gencode arch=compute_20,code=sm_20 -g -G -c GPUHandler.cu CountingHashTable.cu

#### Clean target deletes all generated files ####
clean: 
	rm -rf *.o kmer-counter

# Enable dependency checking
.KEEP_STATE:
.KEEP_STATE_FILE:.make.state.GNU-amd64-Linux

