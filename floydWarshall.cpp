// -*- C++ -*-
//Author : Jae Hong Yoon 
/*
 * floydWarshall.cilk
 * An implementation of blocked Floyed-Warshall algorithm using 
 * Cilk parallelization.
 *
 * Authors: Jaehong Yoon
 * Date: 09 Feb. 2017
 *
 * Copyright (c) belongs to the authors
 *
 */
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
// #include <cilkview.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <functional>
#include <cstring>

#define  __forceinline__ __attribute__((always_inline))

//using namespace std;

/* ==================================================
 *
 * sub2ind - Column-major indexing of 2D arrays
 *
 */
// template <typename T>
// __forceinline__ T sub2ind( T i, T j, T height ) {

  // return (i + height*j);

// }  // end function 'sub2ind'


/* ==================================================
 *
 * initMat - weight matrix (adjacency matrix) initialization
 *
 */
void initMat( double* W, int n) {
			
	cilk_for( int i = 0; i < n; ++i ){
		cilk_for( int j = 0; j < n; ++j ){
			if(i==j){
				W[ i*n + j ] = 0;
			}else{
				int randgen = rand() % 10 + 1;
				if( randgen > 5 ){
					W[ i*n + j ] = 9999999999999999999999999999; // pseudo infinity
				}else{
					W[ i*n + j ] = rand() % 10 + 1;
				}
			}			
		}
	}

	// check the size of array
	std::cout<< "Length of Array = " << (sizeof(W)/sizeof(*W)) << std::endl;
	
	if ( (sizeof(W)/sizeof(*W)) != n*n ){
		std::cerr << "Something went wrong with initialization of matrix"<< std::endl;
		return ;
	}
	
}  // end function 'initMat'


/* ==================================================
 *
 * floydWarshall_calc - return the shortest distance calculated
 * from i to j-th node using k-th node as intermediate node
 *
 */
template<typename T> double floydWarshall_calc(int n, 
											int k, int i, int j,
											T* A){
	// input:
	// n = number of vertex (== ||V||)
	// i,j = node number
	// k = intermediate node number
	// A = adjacency matrix with weight (=W)

	double result = std::min(A[ i*n + j ], 
					A[ i*n + k ]+ A[ k*n + j ]);
	
	return result;
} // end function 'floydWarshall_SP'
 
 
 /* ==================================================
 *
 * floydWarshall_serial - original serial version of Floyed-Warshall
 *
 */
template<typename T> void floydWarshall_serial(int n, T* A){
	// input:
	// n = number of vertex (== ||V||)
	// i,j = node number
	// k = intermediate node number
	// A = adjacency matrix with weight (=W)
	
	int i, j, k;
	for ( k = 0; k < n; ++k ){
		for( i = 0; i < n; ++ i){
			for( j = 0; j < n; ++ j){
				A[ i*n + j ] = floydWarshall_calc(n, k, i, j, A);
			}
		}
	}

} // end function 'floydWarshall_serial'
 
 
 /* ==================================================
 *
 * floydWarshall_SP - simple parallel Floyed-Warshall
 *
 */
template<typename T> void floydWarshall_SP(int n, T* A){
	// input:
	// n = number of vertex (== ||V||)
	// i,j = node number
	// k = intermediate node number
	// A = adjacency matrix with weight (=W)
	
	int i, k;
	for ( k = 0; k < n; ++k ){
		for( i = 0; i < n; ++ i){
			// using cilk_for might work as well since it automatically
			// asks the runtime-API to check if the for loop can be done
			// in parallel. However, it is not guaranteed that it will run
			// parallel. Using cilk_spawn is, therefore, more strict.
			cilk_for( int j = 0; j < n; ++ j){
				A[ i*n + j ] = floydWarshall_calc(n, k, i, j, A);
			}
			// cilk_sync;
		}
	}

} // end function 'floydWarshall_SP'
 
 
/* ==================================================
 *
 * floydWarshall_divide - recursive floydWarshall function which divides
 * A (N X N) to four (N/2 x N/2) matrix until it reaches single element.
 * If the program eaches the element, then floydWarshall_calc is called.
 * In this way, by spawning the threads, we can make the process go parallel.
 *
 * W^k = W^(k/2)*W&(k/2)
 *
 * This code will be used to calculate the all-pair shortest distance of then
 * diagonal submatrix (A_kk)
 */
 
 template<typename T> void floydWarshall_divide(int n,
											int b,
											int k, int i, int j,
											T* A){
												
	// input:
	// n = number of vertex (== ||V||)
	// b = size of a block. This is passed to recursive call
	//     to jump from a block to block for parallel
	//	   computation.
	// i,j = node number
	// k = intermediate node number
	// A = adjacency matrix with weight (=W)											
												
	int bb = b/2; // since b is int type data, this will have a 
				  // flooring effect
	
	if ( b == 1 ){
		A[ i*n + j ] = floydWarshall_calc(n, k, i, j, A);
	}else{
		// divide the adjacency matrix into four submatrix
		
		// calculate APSD accordingly for the A_11 (top & left
		// submatrix of A)
		cilk_spawn floydWarshall_divide(n, bb, k, i, j, A);
		cilk_spawn floydWarshall_divide(n, bb, k, i, j+bb, A);
		cilk_spawn floydWarshall_divide(n, bb, k, i+bb, j, A);
		cilk_spawn floydWarshall_divide(n, bb, k, i+bb, j+bb, A);
		
		// synchronize due to dependency
		cilk_sync;
		
		// calculate APSD accordingly for the A_22 (bottom & right
		// submatrix of A)
		cilk_spawn floydWarshall_divide(n, bb, k+bb, i, j, A);
		cilk_spawn floydWarshall_divide(n, bb, k+bb, i, j+bb, A);
		cilk_spawn floydWarshall_divide(n, bb, k+bb, i+bb, j, A);
		cilk_spawn floydWarshall_divide(n, bb, k+bb, i+bb, j+bb, A);
		
		// synchronize due to dependency
		cilk_sync;
	}												
} // end function 'floydWarshall_divide'


/* ==================================================
 *
 * floydWarshall_propagate - This code will be used to broadcast the block
 * parallel process for pivot row/column matrix. The string input str will be used
 * to determine the direction of propagation.
 *
 */

 template<typename T> void floydWarshall_propagate(int n,
											int b,
											int k, int i, int j,
											char* str,
											T* A){
												
	// input:
	// n = number of vertex (== ||V||)
	// b = size of a block. This is passed to recursive call
	//     to jump from a block to block for parallel
	//	   computation.
	// i,j = node number
	// k = intermediate node number
	// A = adjacency matrix with weight (=W)											
												
    std::cout<< "(n, b, k, i, j) = (" << n << ", " <<b <<", " 
			<<k <<", " << i << ", " << j << " )"  << std::endl;
	std::cout<< "direction = " << str << std::endl;
	
	int bb = b/2; // since b is int type data, this will have a 
				  // flooring effect
	
	if ( b == 1 ){

		A[ i*n + j ] = floydWarshall_calc(n, k, i, j, A);
		std::cout << "calculation for (i, j, k) = ( " << i << ", " << j << ", " 
									<< k << " )" << "done." << std::endl;
	}else{
		std::cout<< "since size of subarray is not 1!!!" << std::endl;
		// case for pivot row matrix
		if(strcmp(str, "column")){
			// the function will be called for pivot row matrix which is not
			// diagonal (A_ij| i~= j). Therefore, it should start propagating on the
			// main direction (column) but also in the other direction (row) to cover
			// the whole submatrix.
			// The order of propagation will be only chaged for row propagation
			// (1) call self recursively
			// (2) broadcast row direction by bb
            std::cout<< "--------------------------------------------" << std::endl;
            std::cout<< "(i, j) = (" << i << ", " << j << ")" << std::endl;
            std::cout<< "(i+bb, j) = (" << i+bb << ", " << j << ")" << std::endl;

			cilk_spawn floydWarshall_propagate(n, bb, k, i, j, "column", A);
			cilk_spawn floydWarshall_propagate(n, bb, k, i+bb, j, "column", A);
			
			cilk_sync;
			
			// (3) now propagate towards column and diagonal
			// due to dependency of the propagated matrix to the pivot column and 
			// non-pivot part of submatrix, these part will be fully calculated
            std::cout<< "(i, j+bb) = (" << i << ", " << j+bb << ")" << std::endl;
            std::cout<< "(i+bb, j+bb) = (" << i+bb << ", " << j+bb << ")" << std::endl;

			cilk_spawn floydWarshall_divide(n, bb, k, i, j+bb, A);
			cilk_spawn floydWarshall_divide(n, bb, k, i+bb, j+bb, A);
			
			cilk_sync;
			
			// (4) now do the calculation for the new intermediate
                        std::cout<< "(i+bb, j) = (" << i+bb << ", " << j << ")" << std::endl;
                        std::cout<< "(i+bb, j+bb) = (" << i+bb << ", " << j+bb << ")" << std::endl;

			cilk_spawn floydWarshall_propagate(n, bb, k+bb, i+bb, j, "column", A);
			cilk_spawn floydWarshall_propagate(n, bb, k+bb, i+bb, j+bb, "column", A);
			
			cilk_sync;
			
			// (5) now calculate the apsd for the original position with new
			// intermediate.
			// Since all dependent matrix are calculate, the values for these can be
			// fully recoverd.
                        std::cout<< "(i, j) = (" << i << ", " << j << ")" << std::endl;
                        std::cout<< "(i, j+bb) = (" << i << ", " << j+bb << ")" << std::endl;

			cilk_spawn floydWarshall_divide(n, bb, k+bb, i, j, A);
			cilk_spawn floydWarshall_divide(n, bb, k+bb, i, j+bb, A);
			
			cilk_sync;
			
		}else{// case for pivot column matrix
			// basic strategy is the same except now the propagation is towards 
			// row direction
            std::cout<< "xxxxx" << std::endl;
            std::cout<< "(i, j) = (" << i << ", " << j << ")" << std::endl;
            std::cout<< "(i+bb, j) = (" << i+bb << ", " << j << ")" << std::endl;

			cilk_spawn floydWarshall_propagate(n, bb, k, i, j, "row", A);
			cilk_spawn floydWarshall_propagate(n, bb, k, i+bb, j, "row", A);
			
			cilk_sync;
			
            std::cout<< "(i, j+bb) = (" << i << ", " << j+bb << ")" << std::endl;
            std::cout<< "(i+bb, j+bb) = (" << i+bb << ", " << j+bb << ")" << std::endl;

			cilk_spawn floydWarshall_divide(n, bb, k, i, j+bb, A);
			cilk_spawn floydWarshall_divide(n, bb, k, i+bb, j+bb, A);
			
			cilk_sync;
			
            std::cout<< "(i, j+bb) = (" << i << ", " << j+bb << ")" << std::endl;
            std::cout<< "(i+bb, j+bb) = (" << i+bb << ", " << j+bb << ")" << std::endl;

			cilk_spawn floydWarshall_propagate(n, bb, k+bb, i, j+bb, "row", A);
			cilk_spawn floydWarshall_propagate(n, bb, k+bb, i+bb, j+bb, "row", A);
			
			cilk_sync;
			
            std::cout<< "(i, j) = (" << i << ", " << j << ")" << std::endl;
            std::cout<< "(i+bb, j) = (" << i+bb << ", " << j << ")" << std::endl;

			cilk_spawn floydWarshall_divide(n, bb, k+bb, i, j, A);
			cilk_spawn floydWarshall_divide(n, bb, k+bb, i+bb, j, A);
			
			cilk_sync;
		}

	}												
} // end function 'floydWarshall_propagate'


/* ==================================================
 *
 * floydWarshall_WU - WU stands for Wrap up. The function pass the top/left
 * submatrix (pivot matrix) first, and then the pivot row/colmn simultaneously.
 * Then, the new pivot matrix is relaxed.
 * 
 * To put it all together, the process is back propagated.
 *
 * ref: Srinivasan, Thanukrishnan, et al. "A scalable parallelization of
 * all-pairs shortest path algorithm for a high performance cluster environment." 
 * Parallel and Distributed Systems, 2007 International Conference on. Vol. 2. 
 * IEEE, 2007.
 *
 */
 
 template<typename T> void floydWarshall_WU(int n,
											int b,
											int k, int i, int j,
											T* A
											){

	// input:
	// n = number of vertex (== ||V||)
	// b = size of a block. This is passed to recursive call
	//     to jump from a block to block for parallel
	//	   computation.
	// i,j = node number
	// k = intermediate node number
	// A = adjacency matrix with weight (=W)
	// str = either 'row' or 'col' to assign direction of broadcasting
												
	int bb = b/2; // since b is int type data, this will have a 
				  // flooring effect
	
	if ( b == 1 ){

//	__cilkrts_set_param("nworkers", "1");

		A[ i*n + j ] = floydWarshall_calc(n, k, i, j, A);
	}else{
	    
		std::cout<< "i = "<< i << ", j= " <<  j << ", k= "<< k << std::endl;
		std::cout<< "n = "<< n << ", b= " <<  b << ", bb= "<< bb << std::endl;
	
		// forward propagation
		// pass pivot matrix (recursively)
		floydWarshall_WU(n, bb, k, i, j, A);			 
												 							 
		// pass pivot row/column. This step can be done paralle
		cilk_spawn floydWarshall_propagate(n, bb, k, i+bb, j, "row", A);
		cilk_spawn floydWarshall_propagate(n, bb, k, i, j+bb, "column", A);
		
		cilk_sync;
		
		// pass the next pivot matrix
		floydWarshall_divide(n, bb, k, i+bb, j+bb, A);
				
		// backward propagation
		// pass next pivot matrix with new intermediate
		floydWarshall_WU(n, bb, k+bb, i+bb, j+bb, A);
				
		// pass the row/column for new intermediate
		cilk_spawn floydWarshall_propagate(n, bb, k+bb, i+bb, j, "column", A);
		cilk_spawn floydWarshall_propagate(n, bb, k+bb, i, j+bb, "row", A);
		
		cilk_sync;
		
		// pass the pivot matrix with new intermediate
		floydWarshall_divide(n, bb, k+bb, i, j, A);
		
		
	}
} // end function 'floydWarshall_WU'


// A simple test harness 
int Floyed_Warshall_Main(int n)
{
	// i, j for for-loop
	int i, j;
	
	// input n = number of vertex in graph (||V||)
	// allocate double type weight matrix (adjacency matrix) with size of
	// (n x n)
	double* W = (double*) calloc(n*n, sizeof(double));
	double* W2 = (double*) calloc(n*n, sizeof(double));
	double* W3 = (double*) calloc(n*n, sizeof(double));
	
//    cilk::cilkview cv; //not sure about this part but uncommented it from
					   // Alex's code
	
	// initialize the weight matrixs
	initMat( W, n );
	W2 = W; W3 = W;


    std::cout << "Calculting all-pair shortest distance for size " << n 
													<< " matrix" << std::endl;

//	cv.start();
    floydWarshall_WU(n, n, 0, 0, 0, W);
//    cv.stop();
//    cv.dump("floydWarshall-results", false);
//	std::cout << "Block version Parallel Floyed-Warshall computation took" << 
//				cv.accumulated_milliseconds() / 1000.f << " seconds" << std::endl;

//	cv.reset();
//	cv.start();
    floydWarshall_SP(n, W2);
//    cv.stop();
//   cv.dump("floydWarshall-results", false);
//	std::cout << "Simple Floyed-Warshall computation took" << 
//				cv.accumulated_milliseconds() / 1000.f << " seconds" << std::endl;
	
//	cv.reset();
//	cv.start();
    floydWarshall_serial(n, W3);
//    cv.stop();
//    cv.dump("floydWarshall-results", false);
//	std::cout << "Serial Floyed-Warshall computation took" << 
//				cv.accumulated_milliseconds() / 1000.f << " seconds" << std::endl;
	
	// confirm that shortest distance calculated for three method matches
	// since the parallel version is the simplest and most accurate, it will be
	// used as reference for validation.
	
	// Block Parallel vs. Serial
	for (i = 0; i < n; ++i){
		for (j = 0; j < n; ++j){
			if(W[ i*n + j ] != W3[ i*n + j ]){
				std::cout<< "Block Parallel Floyed-Warshall failed at (i, j) = (" << i 
					<< ", " << j << ")" << std::endl;
					delete[] W, W2, W3;
					return 1;
			}
		}
	}
	
	// Simple Parallel vs. Serial
	for (i = 0; i < n; ++i){
		for (j = 0; j < n; ++j){
			if(W2[ i*n + j ] != W3[ i*n + j ]){
				std::cout<< "Simple Parallel Floyed-Warshall failed at (i, j) = (" << i 
					<< ", " << j << ")" << std::endl;
					delete[] W, W2, W3;
					return 1;
			}
		}
	}
	
	std::cout << "Floyed-Warshall succeeded for all three implementations." << std::endl;
	delete[] W, W2, W3;
    return 0;
}


int main()
{
	// ask for size of adjacency matrix
	int n;
    std::cout << "Please enter the size of matrix (number of vertex): "<<std::endl;
	std::cin >> n;
	
	if (n <= 0){
		std::cerr << "Invalid argument"<< std::endl;
        std::cerr << "N = number of vetrex in graph" << std::endl;
		std::cerr << "Therefore, N should be greater than 0" << std::endl;
	}

    return Floyed_Warshall_Main(n);
}
