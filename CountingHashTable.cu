#include "CountingHashTable.h"
#include <inttypes.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include "KMerSizes.h"
#include <thread>
#include <chrono>
#include <fstream>

const uint32_t MAX_TRIALS = 60;

template<uint32_t key_size, uint32_t threads_per_block>
__global__ void StartInsertionKernel(KmerKeyValue<key_size>* d_input, uint32_t no_of_keys,
		KmerKeyValue<key_size>* kmer_db, uint64_t kmer_db_max_record_count, uint64_t* cuda_counters) {
	//printf("Executing kernel, about to handle %i keys\n", no_of_keys);

	uint32_t thread_id = (blockIdx.x * blockDim.x + threadIdx.x);
	if (thread_id >= no_of_keys) {
		return;
	}
	//uint64_t index = thread_id * sizeof(KmerKeyValue<key_size>);


	//printf("Printing local entry %"PRIu64", %"PRIu64", %"PRIu64"\n", *(uint64_t*)&d_input[thread_id].getKey(), d_input[thread_id].getCount(), d_input[thread_id].getKey().hash(0));


//    const KmerKeyValue<key_size>* p = &input[index];
//    printf("current thread index:%"PRIu64", key:%"PRIu64"\n", thread_id, *(uint64_t*)&(d_input[thread_id].getKey()));


	const KmerKeyValue<key_size>* key_value = &(d_input[thread_id]);
	bool success = false;
	for(uint32_t trial = 0; trial < MAX_TRIALS && !success; trial++) {

		const KmerKey<key_size>& key = key_value->getKey();
		if (key.isAllA()) {
//			printf("=========================ALL AAAAAAAAAAAAA input_count:%"PRIu64", countbeforeadd=%"PRIu64"\n",
//					d_input[thread_id].getCount(), kmer_db[kmer_db_max_record_count - 1].getCount());
			atomicAdd((unsigned long long int*)(((char*)&kmer_db[kmer_db_max_record_count - 1]) + sizeof(KmerKey<key_size>)), (unsigned long long int)d_input[thread_id].getCount());
//			printf("=========================ALL AAAAAAAAAAAAA count after:%"PRIu64"\n", kmer_db[kmer_db_max_record_count - 1].getCount());
			//atomicAdd(((unsigned long long int*)&cuda_counters[0]), (unsigned long long int)1);
			success = true;
			break;
		}

		uint64_t location = key_value->getKey().hash(trial) % (kmer_db_max_record_count - 1);

//		printf("Generating hash current thread index:%"PRIu64", key:%"PRIu64", hash=%"PRIu64", trial=%i\n",
//				thread_id, *(uint64_t*)&(d_input[thread_id].getKey()), key->hash(trial), trial);

		unsigned long long int zero = 0;
		uint64_t oldValue = atomicCAS((unsigned long long int*)&kmer_db[location], zero, *(unsigned long long int*)&d_input[thread_id].getKey());

//		printf("==========old value:%"PRIu64", new value:%"PRIu64"\n", oldValue, *(unsigned long long int*)&d_input[thread_id].getKey());

		if (oldValue == 0) {
			// a new entry
			success = true;
			for (uint32_t index = sizeof(uint64_t); index < sizeof(KmerKey<key_size>); index++) {
				*(((char*)&kmer_db[location]) + index) = *((char*)&d_input[thread_id].getKey() + index);
			}

			atomicAdd((unsigned long long int*)(((char*)&kmer_db[location]) + sizeof(KmerKey<key_size>)), (unsigned long long int)d_input[thread_id].getCount());
//			atomicAdd(((unsigned long long int*)&kmer_db[location] + sizeof(KmerKey<key_size>)), (unsigned long long int)d_input[thread_id].getCount());
			//unsigned long long int count = 1;
			//atomicAdd((unsigned long long int*)kmer_db_record_count, count);
			//kmer_db_record_count++;
			atomicAdd(((unsigned long long int*)&cuda_counters[0]), (unsigned long long int)1);
//			printf("Inserting new value key:%"PRIu64"\n", *(uint64_t*)&(d_input[thread_id].getKey()));

			break;
		} else if (oldValue == *(uint64_t*)&d_input[thread_id].getKey()) {
			// check whether keys are equal or else continue
			if (kmer_db[location].getKey() == d_input[thread_id].getKey()) {
				atomicAdd((unsigned long long int*)(((char*)&kmer_db[location]) + sizeof(KmerKey<key_size>)), (unsigned long long int)d_input[thread_id].getCount());
				//atomicAdd(((unsigned long long int*)&kmer_db[location] + sizeof(KmerKey<key_size>)), (unsigned long long int)d_input[thread_id].getCount());
				success = true;

//				printf("Key already existing updating key:%"PRIu64"\n", *(uint64_t*)&(d_input[thread_id].getKey()));
			}
		}
	}

	if (!success) {
		// failed to enter has table, hence insert in to failed keys

//		printf("Failed to insert key : %"PRIu64", FialedCount=%"PRIu64"\n", *(uint64_t*)&key, cuda_counters[1]);

		unsigned long long int zero = 0;
		atomicAdd(((unsigned long long int*)&cuda_counters[1]), (unsigned long long int)1);
//		while(atomicCAS((unsigned long long int*)&kmer_db[location], zero, *(unsigned long long int*)&d_input[thread_id].getKey()) != 0) {
//
//		}
	}
}

template<uint32_t key_size>
__device__ inline int32_t firstFreePosition(KmerKeyValue<key_size>* d_input, uint32_t keys_per_thread, int32_t start) {
//	printf("=================firstFreePosition kpt=%" PRIu32 ", start=%" PRIu32 "\n", keys_per_thread, start);
	int32_t i = start;
	for (; i < keys_per_thread; i++) {
		KmerKeyValue<key_size>* kvp = &d_input[i];
		if (kvp->getCount() == 0) {
			break;
		}
	}
	return i;
}

template<uint32_t key_size>
__device__ inline int32_t lastOccupiedPosition(KmerKeyValue<key_size>* d_input, uint32_t keys_per_thread, int32_t last) {
//	printf("=================lastOccupiedPosition kpt=%" PRIu32 ", last=%" PRIu32 "\n", keys_per_thread, last);
	int32_t i = last;
	for (; i > 0; i--) {
		KmerKeyValue<key_size>* kvp = &d_input[i];
		if (kvp->getCount() != 0) {
			break;
		}
	}
	return i;
}

template<uint32_t key_size>
__global__ void ShrinkAndPackHashTable(KmerKeyValue<key_size>* d_input, char* d_output, uint64_t* output_count,
		uint32_t keys_per_thread, uint64_t kmer_db_max_record_count) {
	uint64_t thread_id = (blockIdx.x * blockDim.x + threadIdx.x);

	uint64_t index = thread_id * keys_per_thread;
	if (index > kmer_db_max_record_count) {
		return;
	}

//	printf("=================first index=%" PRIu64 ", last index=%" PRIu64 "\n", index, index + keys_per_thread - 1);

//	printf("==========================thread id =%" PRIu64 ", bloc id=%" PRIu32 ", block dim=%" PRIu32 ", index=%" PRIu32 "\n",
//				thread_id, blockIdx.x, blockDim.x, index);

//	if (thread_id == 0) {
//		printf("=================calling from thread id 0\n");
//	}

	if (index + keys_per_thread > kmer_db_max_record_count) {
		keys_per_thread = (kmer_db_max_record_count % keys_per_thread) + 1;

//		printf("=================The max thread reach here thread_id=%" PRIu64 ", index=%" PRIu64 ", kmer_db_max_record_count=%" PRIu64 ", keys_per_thread=%" PRIu64 "\n",
//						thread_id, index, kmer_db_max_record_count, keys_per_thread);
	}

//	printf("=================================CreateSortedHostData=%" PRIu64 ", keys_per_thread=%" PRIu64 "\n", index, keys_per_thread);

	uint32_t l = firstFreePosition(d_input + index, keys_per_thread, 0);
	uint32_t r = lastOccupiedPosition(d_input + index, keys_per_thread, keys_per_thread - 1);

//	printf("=================================CreateSortedHostData index=%" PRIu64 ",l=%" PRIu32 ", r=%" PRIu32 ", lv=%" PRIu64 ", rv=%" PRIu64 "\n",
//			index, l, r, *(uint64_t*)&d_input[index + l], *(uint64_t*)&d_input[index + r]);
	//&& ((index + l) <= kmer_db_max_record_count) && ((index + r) <= kmer_db_max_record_count)
//	uint64_t count = 0;
//	if (l > r) {
//		count = l;
//	} else {
		while(l < r) {
//			count++;
			memcpy(&d_input[index + l], &d_input[index + r], sizeof(KmerKeyValue<key_size>));
			memset(&d_input[index + r], 0, sizeof(KmerKeyValue<key_size>));
			l = firstFreePosition(d_input + index, keys_per_thread, l + 1);
			r = lastOccupiedPosition(d_input + index, keys_per_thread, r - 1);
		}
//	}

//	printf("=================last valid index = %" PRIu64 "\n", l);

//	if (count > 0) {
//		printf("=================Count = %" PRIu64 ", index=%" PRIu64 "\n", count, index);
		uint64_t oldValue = atomicAdd((unsigned long long int*)output_count, (unsigned long long int)l);
		memcpy(d_output + oldValue * sizeof(KmerKeyValue<key_size>), (char*)(d_input + index), sizeof(KmerKeyValue<key_size>) * l);
//	}

	__syncthreads();
//	printf("=================last valid index = %" PRIu64 ", oldValue=%" PRIu64 "\n", l, oldValue);
}

class KMer32Comparator1 {
public:
	__device__ __host__
	bool operator()(const KMer32 &c1, const KMer32 &c2) {
//		printf("=================================KMer32Comparator1 l=%" PRIu64 ", r=%" PRIu64 "\n", c1.kmer[0], c2.kmer[0]);
		return c1.kmer[0] < c2.kmer[0];
	}
};

class KMer64Comparator1 {
public:
	__device__ __host__
	bool operator()(const KMer64 &c1, const KMer64 &c2) {
//		printf("=================================KMer64Comparator1 l=%" PRIu64 ", r=%" PRIu64 "\n", c1.kmer[0], c2.kmer[0]);
		for (uint32_t i = 0; i < 2; i++) {
			if (c1.kmer[i] < c2.kmer[i])
				return true;
			else if (c1.kmer[i] > c2.kmer[i])
				return false;
		}
		return false;
	}
};

class KMer96Comparator1 {
public:
	__device__ __host__
	bool operator()(const KMer96 &c1, const KMer96 &c2) {
//		printf("=================================KMer96Comparator1 l=%" PRIu64 ", r=%" PRIu64 "\n", c1.kmer[0], c2.kmer[0]);
		for (uint32_t i = 0; i < 3; i++) {
			if (c1.kmer[i] < c2.kmer[i])
				return true;
			else if (c1.kmer[i] > c2.kmer[i])
				return false;
		}
		return false;
	}
};

template<uint32_t key_size>
char* CreateSortedHostData(KmerKeyValue<key_size>* d_input, char* d_output, uint64_t* output_count, uint64_t kmer_db_max_record_count, uint64_t& size) {
	printf("=================================CreateSortedHostData\n");

//	char* temp = new char[kmer_db_max_record_count * 16];
//	CUDA_CHECK_RETURN(cudaMemcpy(temp, d_input, kmer_db_max_record_count * 16, cudaMemcpyDeviceToHost));
//	CUDA_CHECK_RETURN(cudaDeviceSynchronize());
//	printf("%" PRIu64 " %" PRIu64 "\n", *(uint64_t*)(temp + (kmer_db_max_record_count - 1) * 16), *(uint64_t*)(temp + ((kmer_db_max_record_count -1) * 16) + 8));

//	uint32_t threads_per_block = 512;
//	uint32_t block_count = (kmer_db_max_record_count + threads_per_block - 1) / threads_per_block;
//	uint32_t keys_per_thread = 1;
//	if (block_count > 65535) {
//		uint32_t temp = block_count + 65534 / 65535;
//		keys_per_thread = temp;
//		block_count = (kmer_db_max_record_count + keys_per_thread - 1) / keys_per_thread;
//	}

	//=================================================
////	printf("+++++++++++++++++++++++++++++++++++++++++++Before Shrink\n");
//	char* temp1 = new char[kmer_db_max_record_count * sizeof(KmerKeyValue<key_size>)];
//	CUDA_CHECK_RETURN(cudaMemcpy(temp1, d_input, kmer_db_max_record_count * sizeof(KmerKeyValue<key_size>), cudaMemcpyDeviceToHost));
//	CUDA_CHECK_RETURN(cudaDeviceSynchronize());
//	for (uint64_t i = 0; i < kmer_db_max_record_count * sizeof(KmerKeyValue<key_size>); i += sizeof(KmerKeyValue<key_size>)) {
//		if (*(uint64_t*)(temp1 + i + sizeof(KmerKeyValue<key_size>) - sizeof(uint64_t)) > 0) {
//			printf("index=%" PRIu64 ", val=%" PRIu64 ", count=%" PRIu64 "\n", i, *(uint64_t*)(temp1 + i), *(uint64_t*)(temp1 + i + sizeof(KmerKeyValue<key_size>) - sizeof(uint64_t)));
//		}
////		uint64_t v = 1;
////		memcpy(temp1 + i, &i, 8);
////		memcpy(temp1 + i + 8, &v, 8);
//	}
//	CUDA_CHECK_RETURN(cudaMemcpy(d_input, temp1, kmer_db_max_record_count * 16, cudaMemcpyHostToDevice));
//	CUDA_CHECK_RETURN(cudaDeviceSynchronize());
	//=================================================

	uint32_t threads_per_block = 512;// 512;
	uint32_t block_count = 2048;// 2048;
	uint32_t keys_per_thread = kmer_db_max_record_count / (threads_per_block * block_count) + 1;
			// (((kmer_db_max_record_count + block_count - 1) / block_count) + threads_per_block - 1) / threads_per_block;


	//= (((kmer_db_max_record_count + threads_per_block - 1) / threads_per_block) + block_count - 1) / block_count;
			//(kmer_db_max_record_count + block_count - 1) / (block_count + threads_per_block - 1) / threads_per_block;
	//CUDA_CHECK_RETURN(cudaDeviceSynchronize());
	ShrinkAndPackHashTable<<<block_count, threads_per_block>>>(d_input, d_output, output_count, keys_per_thread, kmer_db_max_record_count);
	CUDA_CHECK_RETURN(cudaDeviceSynchronize());

	//=================================================
//	std::this_thread::sleep_for(std::chrono::seconds(20));
//	printf("+++++++++++++++++++++++++++++++++++++++++++After Shrink\n");
//	char* temp = new char[kmer_db_max_record_count * 16];
////	CUDA_CHECK_RETURN(cudaMemcpy(temp, d_output, kmer_db_max_record_count * 16, cudaMemcpyDeviceToHost));
//	CUDA_CHECK_RETURN(cudaMemcpy(temp, d_input, kmer_db_max_record_count * 16, cudaMemcpyDeviceToHost));
//	CUDA_CHECK_RETURN(cudaDeviceSynchronize());
//	printf("%" PRIu64 " %" PRIu64 "\n", *(uint64_t*)(temp + kmer_db_max_record_count * 16), *(uint64_t*)(temp + (kmer_db_max_record_count * 16) + 8));
//	std::ofstream tempfile("/tmp/temp_data_dump.log");
//	for (uint64_t i = 0; i < kmer_db_max_record_count * 16; i += 16) {
//		if (*(uint64_t*)(temp + i + 8) > 0) {
//			printf("%" PRIu64 " %" PRIu64 "\n", *(uint64_t*)(temp + i), *(uint64_t*)(temp + i + 8));
//		}
//		tempfile.write((char*)(uint64_t*)(temp + i), 8);
//		uint32_t v = (uint32_t)(*(uint64_t*)(temp + i + 8));
//		tempfile.write((char*)(uint32_t*)&v, 4);
//	}
//	tempfile.close();
	//=================================================

	uint64_t count = 0;
	CUDA_CHECK_RETURN(cudaMemcpy(&count, (char*)output_count, sizeof(uint64_t), cudaMemcpyDeviceToHost));

	printf("=================================CreateSortedHostData packed record count=%" PRIu64 ", size KMer32=%" PRIu64 ", key-size=%i\n",
			count, sizeof(KMer64), key_size);

	if (key_size == 1) {
		thrust::sort(thrust::cuda::par, (KMer32*)d_output, ((KMer32*)(d_output)) + count, KMer32Comparator1());
	} else if (key_size == 2) {
		thrust::sort(thrust::cuda::par, (KMer64*)d_output, ((KMer64*)(d_output)) + count, KMer64Comparator1());
	} else if (key_size == 3) {
		thrust::sort(thrust::cuda::par, (KMer96*)d_output, ((KMer96*)(d_output)) + count, KMer96Comparator1());
	}

	CUDA_CHECK_RETURN(cudaDeviceSynchronize());

	char* data = new char[count * sizeof(KmerKeyValue<key_size>)];
	CUDA_CHECK_RETURN(cudaMemcpy(data, d_output, count * sizeof(KmerKeyValue<key_size>), cudaMemcpyDeviceToHost));

//	for (uint64_t i = 0; i < count * sizeof(KmerKeyValue<key_size>); i += sizeof(KmerKeyValue<key_size>)) {
//		if (*(uint64_t*)(data + i + sizeof(KmerKeyValue<key_size>) - sizeof(uint64_t)) > 0) {
//			printf("index=%" PRIu64 ", val=%" PRIu64 ", count=%" PRIu64 "\n", i, *(uint64_t*)(data + i), *(uint64_t*)(data + i + sizeof(KmerKeyValue<key_size>) - sizeof(uint64_t)));
//		}
//	}

	size = count * sizeof(KmerKeyValue<key_size>);
	return data;
}

template<uint32_t key_size>
void InsertToHashTable(KmerKeyValue<key_size>* d_input, uint32_t no_of_keys_per_stream, cudaStream_t stream,
		KmerKeyValue<key_size>* kmer_db, uint64_t kmer_db_max_record_count, uint64_t* cuda_counters) {
	printf("Initiating insertion kernel for stream : %"PRIu64", key-size=%i\n", (uint64_t) stream, key_size);


	const uint32_t threads_per_block = 512;
	uint32_t no_of_blocks = (no_of_keys_per_stream + threads_per_block - 1) / threads_per_block;
	StartInsertionKernel<key_size, threads_per_block> <<<no_of_blocks, threads_per_block, 0, stream>>>(d_input, no_of_keys_per_stream,
			kmer_db, kmer_db_max_record_count, cuda_counters);
	CUDA_CHECK_RETURN(cudaStreamSynchronize(stream));
}


#define EXPORT(x) template void InsertToHashTable<x>(KmerKeyValue<x>* input, uint32_t no_of_keys_per_stream, cudaStream_t stream, KmerKeyValue<x>* kmer_db, uint64_t kmer_db_max_record_count, uint64_t* cuda_counters);
EXPORT(1);
EXPORT(2);
EXPORT(3);
EXPORT(4);

#define EXPORT1(x) template char* CreateSortedHostData(KmerKeyValue<x>* d_input, char* d_output, uint64_t* output_count, uint64_t kmer_db_max_record_count, uint64_t& size);
EXPORT1(1);
EXPORT1(2);
EXPORT1(3);
EXPORT1(4);

///**
// * Host function that prepares data array and passes it to the CUDA kernel.
// */
//
///**
// * CUDA kernel function that reverses the order of bits in each element of the array.
// */
//__global__ void multi(void *data) {
//	*(uint32_t*) (data + 8) = 100 + *(uint32_t*) (data + 8);
//	*(uint32_t*) (data + 20) = 100 + *(uint32_t*) (data + 20);
//	*(uint32_t*) (data + 32) = 100 + *(uint32_t*) (data + 32);
//}
//
//#include <iostream>
//#include <vector>
//using namespace std;
//int main(void) {
//
//	cout << "Creating objects" << endl;
//
//	uint32_t input_count = 1 * 1000 * 1000;
//	uint32_t redord_size = 16;
//	uint32_t input_size = redord_size * input_count;
//	uint64_t k = 0;
//	uint64_t v = 1;
//	char* input_data = new char[input_size];
//	for (int64_t i = 0; i < input_size; i += redord_size) {
//		k++;
//		memcpy(input_data + i, &k, sizeof(uint64_t));
//		memcpy(input_data + i + sizeof(uint64_t), &v, sizeof(uint64_t));
//	}
////	uint64_t x = 16;
////	memcpy(input_data + 32, &x, sizeof(uint64_t));
//
//	cout << "Objection completed" << endl;
//
////	for (int32_t i = 0; i < input_size; i += 12) {
////		KmerKeyValue<1>* p = (KmerKeyValue<1>*) &input_data[i];
////		cout << "host:" << *(uint64_t*) &p->getKey() << "|" << p->getCount() << endl;
////	}
//
//	char* d_data = NULL;
//	CUDA_CHECK_RETURN(cudaMalloc((void** ) &d_data, input_size));
//	CUDA_CHECK_RETURN(cudaMemset(d_data, 0, input_size));
//	CUDA_CHECK_RETURN(cudaMemcpy(d_data, input_data, input_size, cudaMemcpyHostToDevice));
//
////	multi<<<1, 1>>>(d_data);
////
////	CUDA_CHECK_RETURN(cudaMemcpy(input_data, d_data, 32, cudaMemcpyDeviceToHost));
////
////	for (int i = 0; i < 36; i++) {
////		cout << "Dev:" << (int32_t)input_data[i] << endl;
////	}
//
//
//	cout << "Invoking HT" << endl;
//	CountingHashTable<1>* cht = new CountingHashTable<1>(0, 1, 50 * 1000 * 1000, 1000);
//
//	bool capable = false;
//	capable = cht->Insert((KmerKeyValue<1>*) d_data, input_count);
//	cout << "HT Completed" << "|" << capable << endl;
//
//	for (int i = 0; i < 5; i++) {
//		capable = cht->Insert((KmerKeyValue<1>*) d_data, input_count);
//	    cout << "HT Completed" << "|" << capable << "|" << i << endl;
//	}
//
////	const int32_t xx = 1;
////	vector<KmerKeyValue<1>> table_data;
////	cht->Extract(table_data);
////	cout << table_data.size() << endl;
//
//	cout << "TB size : " << cht->GetSuccessCount() << "|" << cht->GetFailureCount() << endl;
//
//	cht->Dump();
//
//	return 0;
//}

