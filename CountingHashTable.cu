#include "CountingHashTable.h"
#include <inttypes.h>

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

		if (key_value->getKey().isAllA()) {
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
			for (uint32_t index = sizeof(uint64_t); index < sizeof(KmerKey<key_size>); index += sizeof(uint64_t)) {
				*(((uint64_t*)&kmer_db[location]) + index) = *((uint64_t*)&d_input[thread_id].getKey() + index);
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
void InsertToHashTable(KmerKeyValue<key_size>* d_input, uint32_t no_of_keys_per_stream, cudaStream_t stream,
		KmerKeyValue<key_size>* kmer_db, uint64_t kmer_db_max_record_count, uint64_t* cuda_counters) {
	printf("Initiating insertion kernel for stream : %"PRIu64"\n", (uint64_t) stream);


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

