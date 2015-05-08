/*****************************************************************************/
// nvcc -O1 -o bpsw bpsw.cu -lrt -lm

#include <cstdio>
#include <cstdlib>
#include <math.h>

// Assertion to check for errors
#define CUDA_SAFE_CALL(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
	if (code != cudaSuccess) 
	{
		fprintf(stderr,"CUDA_SAFE_CALL: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

#define TILE_WIDTH 						128
#define NUM_ELEMENTS					10240						//change this to change how many elements to process
#define START_VAL						170101
#define PRINT_TIME 						1



void initializeArray1D(long *arr, long len, long start); 
void initializeArrayConst(int *arr, long len);
void initializeArraySPrimes(long *arr);
int pruneArray(long *arrSrc, long *arrDest, int *result, long len);

__constant__ long d_sPrimes[256];							//declare constant array for small primes

__global__ void kernel_trialDiv (long* n, int* r) {
	int bx = blockIdx.x;      // ID thread
	int tx = threadIdx.x; 
	int i=0;

	// Identify the row and column of the Pd element to work on
	long memIndex = bx*TILE_WIDTH+tx;
	for (i = 0; i < 256; i++)
	{
//		r[memIndex] = ((n[memIndex])%(d_sPrimes[i]) == 0)? (r[memIndex] - 1) : r[memIndex];			//ternary is slower than if statement
		if (n[memIndex] % d_sPrimes[i] == 0)
			r[memIndex] = r[memIndex] - 1;															//r decreases from 1. Only 1s are prime candidates
	}

	__syncthreads();
}

__global__ void kernel_jacobi(long* nArray, long* dArray, long len) {
	int bx = blockIdx.x;      // ID thread
	int tx = threadIdx.x;
	int result, t;
	long d, dAbs, sign, temp, n1, d1;
	// Identify the row and column of the Pd element to work on
	long memIndex = bx*TILE_WIDTH + tx;
	if (memIndex < len)							//out of bounds checking - some threads will be doing nothing
	{
		result = 0;
		dAbs = 5;
		sign = 1;

		while (result != -1)				//if result != -1, increment d and try again
		{
			n1 = nArray[memIndex];				//reinitialize n1 to n
			d = dAbs*sign;
			t = 1;
			d1 = d;							//reinitialize d1 to d
			d1 = d1 % n1;

			while (d1 != 0)
			{
				while (d1 % 2 == 0)        //while d is even 
				{
					d1 = d1 / 2;
					if (n1 % 8 == 3 || n1 % 8 == 5) t = -t;
				}
				temp = d1;
				d1 = n1;
				n1 = temp;
				if ((d1 % 4 == 3) && (n1 % 4 == 3)) t = -t;
				d1 = d1 % n1;
			}
			if (n1 == 1) result = t;
			else result = 0;
			dAbs = dAbs + 2;
			sign = sign * -1;
		}
	}
	__syncthreads();
	if (memIndex < len)
		dArray[memIndex] = d;
	__syncthreads();
}

__global__ void kernel_lucas(long* nArray, long* dArray, int* rArray, long len) {
	int bx = blockIdx.x;      // ID thread
	int tx = threadIdx.x;
	int i, length;
	long long d, n;
	long long q, q2, u, u2, uold, v, v2, t;

	// Identify the row and column of the Pd element to work on
	long memIndex = bx*TILE_WIDTH + tx;
	if (memIndex < len)							//out of bounds checking - some threads will be doing nothing
	{
		d = (long long) dArray[memIndex];
		n = (long long) nArray[memIndex];
		q = (1 - d) / 4;
		u = 0;
		v = 2;
		u2 = 1;
		v2 = 1;
		q2 = 2 * q;
		t = (n + 1) / 2;						//theta
		length = 32 - __clz(t); //length of our number in bits. //clz(b00010010) = 3 	

		for (i = 0; i < length; i++)
		{
			u2 = (u2 * v2) % n;
			v2 = (v2 * v2 - q2) % n;
			if (t & 1)				//mask = 1
			{
				uold = u;
				u = (u2 * v) + (u * v2);
				u = (u % 2 == 1) ? u + n : u;
				u = (u / 2) % n;
				v = (v2 * v) + (u2 * uold * d);
				v = (v % 2 == 1) ? v + n : v;
				v = (v / 2) % n;
			}

			q = (q*q) % n;
			q2 = q + q;

			t = t >> 1;
		}
		
	}
	__syncthreads();
	if (memIndex < len)
		rArray[memIndex] = (u == 0);

}

int main(int argc, char **argv){
	int i, j;
	long arrLen = 0;
	long start = START_VAL;
		
	// GPU Timing variables
	cudaEvent_t start_program, stop_program;
	cudaEvent_t start_tDiv, stop_tDiv; 
	cudaEvent_t start_memOp, stop_memOp;
	cudaEvent_t start_jac, stop_jac;
	cudaEvent_t start_luc, stop_luc;
	


	float elapsed_gpu;
	
	// Arrays on GPU global memory
	long *d_n1, *d_n3, *d_d;
	int *d_r1, *d_r3;

	// Arrays on the host memory
	long *h_n1, *h_n3, *h_n4;  
	int *h_r1, *h_r3;
//	long *h_sPrimes;	
	
	if (argc > 1) {
		arrLen  = atoi(argv[1]);
	}
	else {
		arrLen = NUM_ELEMENTS;											//arrLen = total number of elements to process
	}
	
	long gridSize = arrLen / TILE_WIDTH;								//gridSize = number of blocks to use
	printf("Number of numbers to check = %d\n", arrLen);

	// Allocate GPU memory
	long allocSizeL = arrLen * sizeof(long);
	long allocSizeInt = arrLen * sizeof(int);
	long allocSizeSP = 256 * sizeof(long);
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_n1, allocSizeL));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_r1, allocSizeInt));

	
	// Allocate arrays on host memory
	h_n1 = (long *) malloc(allocSizeL);
	h_r1 = (int *) malloc(allocSizeInt);
//	h_sPrimes = (long *) malloc(allocSizeSP);
	

	
	// Initialize the host arrays
	printf("\nInitializing the arrays ...");
	initializeArray1D(h_n1, arrLen, start);
	initializeArrayConst(h_r1, arrLen);
	long h_sPrimes[256] =										//First 256 primes
	{						
		2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
		31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
		73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
		127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
		179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
		233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
		283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
		353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
		419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
		467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
		547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
		607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
		661, 673, 677, 683, 691, 701, 709, 719, 727, 733,
		739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
		811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
		877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
		947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013,
		1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
		1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
		1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223,
		1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291,
		1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373,
		1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451,
		1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511,
		1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583,
		1597, 1601, 1607, 1609, 1613, 1619
	};

	printf("\t... done\n\n");

#if PRINT_TIME
	// Create the cuda events for mem ops
	cudaEventCreate(&start_program);
	cudaEventCreate(&stop_program);
	// Record event on the default stream
	cudaEventRecord(start_program, 0);
#endif


#if PRINT_TIME
	// Create the cuda events for mem ops
	cudaEventCreate(&start_memOp);
	cudaEventCreate(&stop_memOp);
	// Record event on the default stream
	cudaEventRecord(start_memOp, 0);
#endif

	// declare shape of block - 256x1 block, 64x1 grid
	dim3 dimBlock1(TILE_WIDTH,1);
	dim3 dimGrid1(gridSize,1);

	// Transfer the arrays to the GPU memory
	CUDA_SAFE_CALL(cudaMemcpy(d_n1, h_n1, allocSizeL, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_r1, h_r1, allocSizeInt, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_sPrimes, h_sPrimes, allocSizeSP));	


	printf("\nPotential primes left = %ld\n", arrLen);

#if PRINT_TIME
	// Create the cuda events for trial division
	cudaEventCreate(&start_tDiv);
	cudaEventCreate(&stop_tDiv);
	// Record event on the default stream
	cudaEventRecord(start_tDiv, 0);
#endif

  	
	// Launch the kernel
	kernel_trialDiv<<<dimGrid1, dimBlock1>>>(d_n1, d_r1);

	
#if PRINT_TIME
	// Stop and destroy the timer for trial division
	cudaEventRecord(stop_tDiv,0);
	cudaEventSynchronize(stop_tDiv);
	cudaEventElapsedTime(&elapsed_gpu, start_tDiv, stop_tDiv);
	printf("\nGPU time - trial division kernel only: %f (msec)\n", elapsed_gpu);
	cudaEventDestroy(start_tDiv);
	cudaEventDestroy(stop_tDiv);
#endif
	
	
	// Check for errors during launch
	CUDA_SAFE_CALL(cudaPeekAtLastError());
	
	// Transfer the results back to the host
	CUDA_SAFE_CALL(cudaMemcpy(h_n1, d_n1, allocSizeL, cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(h_r1, d_r1, allocSizeInt, cudaMemcpyDeviceToHost));
	

	
#if PRINT_TIME
	// Stop and destroy the timer for memory operations
	cudaEventRecord(stop_memOp,0);
	cudaEventSynchronize(stop_memOp);
	cudaEventElapsedTime(&elapsed_gpu, start_memOp, stop_memOp);
	printf("\nGPU time - including memory operations: %f (msec)\n", elapsed_gpu);
	cudaEventDestroy(start_memOp);
	cudaEventDestroy(stop_memOp);
#endif

	// Free-up device memory
	CUDA_SAFE_CALL(cudaFree(d_n1));
	CUDA_SAFE_CALL(cudaFree(d_r1));


	// Allocate memory for n3,r3 on host
	h_n3 = (long *) malloc(allocSizeL);
	h_r3 = (int *) malloc(allocSizeInt);

	arrLen = pruneArray(h_n1,h_n3,h_r1,arrLen);						//copy h_n1 to h_n3, only potential primes. arrLen gets length of new array
	printf("\nPotential primes left = %ld\n", arrLen);

	// Free-up host memory
	free(h_n1);
	free(h_r1);

// done with trial division	
////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// calculate jacobi value for each n


	gridSize = (arrLen / TILE_WIDTH) + (arrLen % TILE_WIDTH != 0);				//calculate new number of blocks. Might not be divisible by tile_width, so add 1

	// declare shape of block - 1x64 block, 64x64 grid
	dim3 dimBlock3(TILE_WIDTH, 1);
	dim3 dimGrid3(gridSize, 1);

#if PRINT_TIME
	// Create the cuda events for mem ops
	cudaEventCreate(&start_memOp);
	cudaEventCreate(&stop_memOp);
	// Record event on the default stream
	cudaEventRecord(start_memOp, 0);
#endif


	CUDA_SAFE_CALL(cudaMalloc((void **)&d_n3, allocSizeL));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_d, allocSizeL));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_r3, allocSizeInt));

	CUDA_SAFE_CALL(cudaMemcpy(d_n3, h_n3, allocSizeL, cudaMemcpyHostToDevice));

#if PRINT_TIME
	// Create the cuda events for jacobi
	cudaEventCreate(&start_jac);
	cudaEventCreate(&stop_jac);
	// Record event on the default stream
	cudaEventRecord(start_jac, 0);
#endif


	kernel_jacobi << <dimGrid3, dimBlock3 >> >(d_n3, d_d, arrLen);

#if PRINT_TIME
	// Stop and destroy the timer for jacobi
	cudaEventRecord(stop_jac, 0);
	cudaEventSynchronize(stop_jac);
	cudaEventElapsedTime(&elapsed_gpu, start_jac, stop_jac);
	printf("\nGPU time - jacobi kernel only: %f (msec)\n", elapsed_gpu);
	cudaEventDestroy(start_jac);
	cudaEventDestroy(stop_jac);
#endif

#if PRINT_TIME
	// Stop and destroy the timer for memory operations
	cudaEventRecord(stop_memOp, 0);
	cudaEventSynchronize(stop_memOp);
	cudaEventElapsedTime(&elapsed_gpu, start_memOp, stop_memOp);
	printf("\nGPU time - including memory operations: %f (msec)\n", elapsed_gpu);
	cudaEventDestroy(start_memOp);
	cudaEventDestroy(stop_memOp);
#endif

	// Check for errors during launch
	CUDA_SAFE_CALL(cudaPeekAtLastError());

//done with calculating jacobi value
////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
//run lucas probable prime test

#if PRINT_TIME
	// Create the cuda events for mem ops
	cudaEventCreate(&start_memOp);
	cudaEventCreate(&stop_memOp);
	// Record event on the default stream
	cudaEventRecord(start_memOp, 0);
#endif

#if PRINT_TIME
	// Create the cuda events for lucas
	cudaEventCreate(&start_luc);
	cudaEventCreate(&stop_luc);
	// Record event on the default stream
	cudaEventRecord(start_luc, 0);
#endif

	kernel_lucas << <dimGrid3, dimBlock3 >> >(d_n3, d_d, d_r3, arrLen);

#if PRINT_TIME
	// Stop and destroy the timer for lucas
	cudaEventRecord(stop_luc, 0);
	cudaEventSynchronize(stop_luc);
	cudaEventElapsedTime(&elapsed_gpu, start_luc, stop_luc);
	printf("\nGPU time - lucas kernel only: %f (msec)\n", elapsed_gpu);
	cudaEventDestroy(start_luc);
	cudaEventDestroy(stop_luc);
#endif

	// Check for errors during launch
	CUDA_SAFE_CALL(cudaPeekAtLastError());
	CUDA_SAFE_CALL(cudaMemcpy(h_r3, d_r3, allocSizeInt, cudaMemcpyDeviceToHost));

#if PRINT_TIME
	// Stop and destroy the timer for memory operations
	cudaEventRecord(stop_memOp, 0);
	cudaEventSynchronize(stop_memOp);
	cudaEventElapsedTime(&elapsed_gpu, start_memOp, stop_memOp);
	printf("\nGPU time - including memory operations: %f (msec)\n", elapsed_gpu);
	cudaEventDestroy(start_memOp);
	cudaEventDestroy(stop_memOp);
#endif

#if PRINT_TIME
	// Stop and destroy the timer for program
	cudaEventRecord(stop_program, 0);
	cudaEventSynchronize(stop_program);
	cudaEventElapsedTime(&elapsed_gpu, start_program, stop_program);
	printf("\nGPU time - total time: %f (msec)\n", elapsed_gpu);
	cudaEventDestroy(start_program);
	cudaEventDestroy(stop_program);
#endif

	// Allocate memory for n4 on host
	allocSizeL = arrLen * sizeof(long);
	h_n4 = (long *)malloc(allocSizeL);

	arrLen = pruneArray(h_n3, h_n4, h_r3, arrLen);						//copy h_n3 to h_n4, only primes. arrLen gets length of new array
	printf("\nNumber of prime numbers = %ld\n", arrLen);
	printf("\nFirst 100 primes detected: \n");
	for (i = 0; i < 10; i++)
	{
		for (j = 0; j < 10; j++)
			printf("%ld ", h_n4[10 * i + j]);
		printf("\n");
	}

	free(h_n3);
	free(h_r3);
	CUDA_SAFE_CALL(cudaFree(d_n3));
	CUDA_SAFE_CALL(cudaFree(d_r3));	
	CUDA_SAFE_CALL(cudaFree(d_d));
		
	return 0;

	
}

void initializeArray1D(long *arr, long len, long start) 		//each element increments by 2
{
	long i;
	for (i = 0; i < len; i++)
	{
		arr[i] = start;
		start += 2;
	}
}

void initializeArrayConst(int *arr, long len)				//initialize result matrix to 2
{
	long i;
	for (i = 0; i < len; i++)
		arr[i] = 1;
}

int pruneArray(long *arrSrc, long *arrDest, int *result, long len)
{
	long i;
	long index = 0;					//length of new array
	for (i = 0; i < len; i++)
	{
		if (result[i] == 1)							//only keep potential primes
		{
			arrDest[index] = arrSrc[i];			//write to the new array with current value
			index = index + 1;
		}	
	}
	return (index);
}