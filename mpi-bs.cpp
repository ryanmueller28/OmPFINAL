#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <random>
#include <limits>
#include "mpi.h"

/**
 * Parallelized bucket sort program using MPI.
 *
 * Fall 2019 IST 246 Final Assignment
 */

/**
 * Determine if an array of unsigned integers is sorted.
 *
 * \param pArr Pointer to the first element of the array.
 *
 * \param n Size of the array.
 *
 * \return True if the array is sorted in non-descending order, false
 * otherwise.
 */
bool isSorted(unsigned *pArr, int n) {
	for (int i = 0; i < n - 1; i++) {
		if (pArr[i + 1] < pArr[i]) {
			return false;
		} // if
	} // for
	return true;
}

/**
 * Application entry point.
 *
 * \param argc Number of command line arguments; ignored by this app.
 * 
 * \param pArgv Array of command line arguments; ignored by this app. 
 */
int main(int argc, char* pArgv[]) {
	using namespace std;

	int rank;		// processor's rank ID
	int nProcs;		// number of processors in use

	// initialize MPI constructs
	MPI_Init(&argc, &pArgv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	unsigned *pMaster = 0;
	int masterN = 100000000;

	// root process fills master array
	if (rank == 0) {
		pMaster = new unsigned[masterN];

		mt19937 prng(chrono::system_clock::now().time_since_epoch().count());
		for (int i = 0; i < masterN; i++) {
			pMaster[i] = (unsigned)prng();
		}
	}


	/****************************************************************************
	 * Put your bucketsort code here!                                           *
	 ****************************************************************************/
	/*
	* Bucket Sort by Ryan Mueller, Kira Reisdorph, Jacob Williams, and Devan Standley
	* 
	* 
	* 
	*/
	
	unsigned MSG_NUMBERS = 0;
	//Ryan: define size variables

	int offset = 0;

	//2m + 1 <= nProcs - 1
	//m = tM / 2
	int tM = nProcs - 1;
	int m = tM / 2; // Processors being used at any one time

	/*
	* Ryan:
	* Below variables are for indexing the array
	* by bucket values
	*/
	unsigned prevIndex = ((numeric_limits<unsigned>::max() / m) * (rank - 1));
	unsigned curIndex = ((numeric_limits<unsigned>::max() / m) * rank);
	unsigned distIndex = curIndex - prevIndex;

	int bucketSize = masterN / tM; //MassterN / 2m

	unsigned *pSend = new unsigned[bucketSize];


	unsigned *pRecv = new unsigned[bucketSize];

	//Ryan: first, segment array into segments based on rank
	//Ryan: if at master rank, segment array
	if (rank == 0){
		
		for (int i = 0; i < nProcs; i++){
			MPI_Send(pMaster + offset, bucketSize, MPI::UNSIGNED, i, MSG_NUMBERS, MPI_COMM_WORLD);
			offset += bucketSize; // moves offset value to next bucket size
		}
		
	}else if (rank != 0){
		MPI_Status status;
		offset = 0;
		for (int i = 1; i <= nProcs; i++){
			//error occured MPI_Recv, on communicator: MPI_COMMWORLD, MPI_ERR_BUFFER:Invalid pointer
			MPI_Recv(pMaster + offset, bucketSize, MPI::UNSIGNED, i, MSG_NUMBERS, MPI_COMM_WORLD, &status);
			offset += bucketSize;
		}
		
	}
	
	if (rank % 2 == 0 && rank != 0) {

		//If an even processor: send to odd processors first
		MPI_Status status;
		for (int i = 1; i <= m; i += 2)
		{
			if (i != rank)
			{
				//Ryan: need to look at the tag before COMM_WORLD, unsure if right
				MPI_Send(pSend, bucketSize, MPI::UNSIGNED, i, MSG_NUMBERS, MPI_COMM_WORLD);
			}
		}

		//receive from odd processors
		
		for (int i = 1; i <= m; i += 2){
			if (i != rank){
				//Ryan: need to look at the tag before COMM_WORLD, unsure if right
				MPI_Recv(pRecv, bucketSize, MPI::UNSIGNED, i, MSG_NUMBERS, MPI_COMM_WORLD, &status);
			}
		}
	} else {

		//Ryan: receive segmented array to even processors
		MPI_Status status;
		for (int i = 0; i <= m; i += 2){
			if (i != rank){
				//Ryan: need to look at the tag before COMM_WORLD, unsure if right
				MPI_Recv(pRecv, bucketSize, MPI::UNSIGNED, i, MSG_NUMBERS, MPI_COMM_WORLD, &status);
			}
		}

		//Ryan: send to even processors

		for (int i = 0; i <= m; i += 2){
			if (i != rank){
				//Ryan: need to look at the tag before COMM_WORLD, unsure if right
				MPI_Send(pSend, bucketSize, MPI::UNSIGNED, i, MSG_NUMBERS, MPI_COMM_WORLD);
			}
		}
	}

	/*
	* Ryan:
	* Send data to buckets based on index numbers
	* This portion divides send and receive based on odd and even
	* But sends the numbers to buckets based on index and rank
	* 
	*/


	if(rank == 0) {
		// master validates sorting and frees memory
		printf("**********************************************************\n");
		if (isSorted(pMaster, masterN)) {
			printf("Array is sorted.\n");
		} else {
			printf("Array is not sorted.\n");
		}
		printf("**********************************************************\n");

		delete[] pMaster;

	} 

	MPI_Finalize();

	return EXIT_SUCCESS;
}
