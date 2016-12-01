/*
 * CountingHashTable.h
 *
 *  Created on: Nov 30, 2016
 *      Author: jayanga
 */

#ifndef COUNTINGHASHTABLE_H_
#define COUNTINGHASHTABLE_H_

#include "CountingHashTableBase.h"

template<uint32_t key_size>
class CountingHashTable : public CountingHashTableBase<key_size> {
public:
	CountingHashTable(const uint32_t device_id = 0, const uint32_t no_of_streams = 1,
			const uint32_t kmer_db_size = 1 * 1000 * 1000, const uint32_t kmer_failed_db_size = 1000) :
		CountingHashTableBase<key_size>(device_id, no_of_streams, kmer_db_size, kmer_failed_db_size) {

	}
};

#endif /* COUNTINGHASHTABLE_H_ */
