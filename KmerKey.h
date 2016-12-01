/*
 * KmerKey.h
 *
 *  Created on: Nov 29, 2016
 *      Author: jayanga
 */

#ifndef KMERKEY_H_
#define KMERKEY_H_

#include <stdint.h>

template<uint32_t key_size>
class KmerKey {
private:
	uint64_t _key[key_size];

public:
	__host__ __device__ uint64_t hash(uint32_t trial) const {
//		uint64_t hash = 0;
//		for (uint32_t i = 0; i < key_size; i++) {
//			hash += (_key[i] >> 32) * trial;
//			hash += ((_key[i] << 32) >> 32);
//		}
//		return hash;
		uint8_t* keyBytes = (uint8_t*) &_key;

		uint32_t h1 = 0;
		uint32_t h2 = 0;
		for (uint32_t i = 0; i < key_size * 8; i++) {
			h1 = 31 * h1 + keyBytes[i];
			h2 = 37 * h2 + ~keyBytes[i];
		}
		return h1 + trial * h2;
	}

	__host__ __device__ bool operator<(const KmerKey<key_size>& other) const {
		for (uint32_t i = 0; i < key_size; i++) {
			if (_key[i] < other._key[i])
				return true;
			else if (_key[i] > other._key[i])
				return false;
		}
		return false;
	}

	__device__ __host__ bool operator==(const KmerKey<key_size>& other) const {
		for (uint32_t i = 0; i < key_size; i++) {
			if (_key[i] != other._key[i])
				return false;
		}
		return true;
	}

	__device__ __host__ bool isAllA() const {
		uint64_t temp = 0;
		for (uint32_t i = 0; i < key_size; i++) {
			temp += _key[i];
		}
		return temp == 0;
	}
};

#endif /* KMERKEY_H_ */
