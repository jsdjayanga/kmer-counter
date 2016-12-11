/*
 * KMerCounterUtils.h
 *
 *  Created on: Nov 11, 2016
 *      Author: jayanga
 */

#ifndef KMERCOUNTERUTILS_H
#define KMERCOUNTERUTILS_H

#ifdef DEBUG_BUILD
#define DEBUG(k, v) do { printf("%s %"PRIu64"\n", k, v); } while (0)
#else
#define DEBUG(k, v)
#endif


inline bool lessThan(char* first, char* second, uint32_t kmer_length) {
	if (kmer_length <= 32) {
		return *(uint64_t*) (first) < *(uint64_t*) (second);
	}
	return false;
}

inline bool equals(char* first, char* second, uint32_t kmer_length) {
	if (kmer_length <= 32) {
		return *(uint64_t*) (first) == *(uint64_t*) (second);
	}
	return false;
}

#endif /* KMERCOUNTERUTILS_H */
