/*
 * KMerSizes.h
 *
 *  Created on: Nov 10, 2016
 *      Author: jayanga
 */

#pragma once

struct KMer32 {
	uint64_t kmer[1];
	uint64_t count;
};

struct KMer64 {
	uint64_t kmer[2];
	uint32_t count;
};

struct KMer96 {
	uint64_t kmer[3];
	uint32_t count;
};

struct KMer128 {
	uint64_t kmer[4];
	uint32_t count;
};
