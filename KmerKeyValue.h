/*
 * KmerKeyValue.h
 *
 *  Created on: Nov 29, 2016
 *      Author: jayanga
 */

#ifndef KMERKEYVALUE_H_
#define KMERKEYVALUE_H_

#include <stdint.h>
#include "KmerKey.h"

template<uint32_t key_size>
class KmerKeyValue {
	KmerKey<key_size> _key;
	uint64_t _count; // making count uint64 to get proper memory alignment

public:
	__device__ __host__ KmerKeyValue() :
			_count(0) {
	}

	KmerKeyValue(const KmerKey<key_size>& key, const uint32_t count) :
			_key(key), _count(count) {
	}

	__device__ __host__ uint64_t getCount() const {
		return _count;
	}

	__device__ __host__ uint64_t setCount(uint64_t count) {
			_count = count;
	}

	__device__ __host__ const KmerKey<key_size>& getKey() const {
		return _key;
	}
};

#endif /* KMERKEYVALUE_H_ */
