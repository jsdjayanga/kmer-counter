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


inline bool lessThan(char* first, char* second, uint32_t key_size_in_longs) {
//	if (kmer_length <= 32) {
//		return *(uint64_t*) (first) < *(uint64_t*) (second);
//	}
//	return false;

	for (uint32_t i = 0; i < key_size_in_longs; i++) {
		if (*(uint64_t*)(first + (i * sizeof(uint64_t))) < *(uint64_t*)(second + (i * sizeof(uint64_t))))
			return true;
		else if (*(uint64_t*)(first + (i * sizeof(uint64_t))) > *(uint64_t*)(second + (i * sizeof(uint64_t))))
			return false;
	}
	return false;
}

inline bool equals(char* first, char* second, uint32_t key_size_in_longs) {
//	if (kmer_length <= 32) {
//		return *(uint64_t*) (first) == *(uint64_t*) (second);
//	}
//	return false;

	for (uint32_t i = 0; i < key_size_in_longs; i++) {
		if (*(uint64_t*)(first + (i * sizeof(uint64_t))) != *(uint64_t*)(second + (i * sizeof(uint64_t))))
			return false;
	}
	return true;
}

#endif /* KMERCOUNTERUTILS_H */
