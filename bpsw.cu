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

#define TILE_WIDTH 						64
#define NUM_ELEMENTS					10240						//change this to change how many elements to process
#define START_VAL						170101
#define PRINT_TIME 						1



void initializeArray1D(long long int *arr, int len, long long int start); 
void initializeArrayConst(int *arr, int len);
void initializeArraySPrimes(long long int *arr);
int pruneArray(long long int *arrSrc, long long int *arrDest, int *result, int len);

__constant__ long long int d_sPrimes[256];							//declare constant array for small primes

__global__ void kernel_trialDiv (long long int* n, int* r) {
	int bx = blockIdx.x;      // ID thread
	int tx = threadIdx.x; 
	int i;
	
	// Identify the row and column of the Pd element to work on
	int memIndex = bx*TILE_WIDTH+tx;
	for (i = 0; i < 256; i++)
	{
//		r[memIndex] = ((n[memIndex])%(d_sPrimes[i]) == 0)? (r[memIndex] - 1) : r[memIndex];			//ternary is slower than if statement
		if (n[memIndex] % d_sPrimes[i] == 0)
			r[memIndex] = r[memIndex] - 1;
	}
//	__syncthreads();
}

int main(int argc, char **argv){
	int arrLen = 0;
//	int i;
	long long int start = START_VAL;
		
	// GPU Timing variables
	cudaEvent_t start_tDiv, stop_tDiv; 
	cudaEvent_t start_memOp, stop_memOp;
	float elapsed_gpu;
	
	// Arrays on GPU global memory
	long long int *d_n1, *d_n2;
	int *d_r1, *d_r2;

	// Arrays on the host memory
	long long int *h_n1, *h_n2;  
	int *h_r1, *h_r2;
//	long long int *h_sPrimes;	
	
	if (argc > 1) {
		arrLen  = atoi(argv[1]);
	}
	else {
		arrLen = NUM_ELEMENTS;										//arrLen = total number of elements to process
	}
	
	int gridSize = arrLen / TILE_WIDTH;								//gridSize = number of blocks to use
	printf("Length of the array = %d\n", arrLen);

	// Allocate GPU memory
	int allocSizeLL = arrLen * sizeof(long long int);
	int allocSizeInt = arrLen * sizeof(int);
	int allocSizeSP = 256 * sizeof(long long int);
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_n1, allocSizeLL));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_r1, allocSizeInt));

	
	// Allocate arrays on host memory
	h_n1 = (long long int *) malloc(allocSizeLL);
	h_r1 = (int *) malloc(allocSizeInt);
//	h_sPrimes = (long long int *) malloc(allocSizeSP);
	

	
	// Initialize the host arrays
	printf("\nInitializing the arrays ...");
	initializeArray1D(h_n1, arrLen, start);
	initializeArrayConst(h_r1, arrLen);
	long long int h_sPrimes[256] =										//First 256 primes
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
	cudaEventCreate(&start_memOp);
	cudaEventCreate(&stop_memOp);
	// Record event on the default stream
	cudaEventRecord(start_memOp, 0);
#endif

	// declare shape of block - 1x64 block, 64x64 grid
	dim3 dimBlock(TILE_WIDTH,1);
	dim3 dimGrid(gridSize,1);

	// Transfer the arrays to the GPU memory
	CUDA_SAFE_CALL(cudaMemcpy(d_n1, h_n1, allocSizeLL, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_r1, h_r1, allocSizeInt, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_sPrimes, h_sPrimes, allocSizeSP));	


#if PRINT_TIME
	// Create the cuda events for trial division
	cudaEventCreate(&start_tDiv);
	cudaEventCreate(&stop_tDiv);
	// Record event on the default stream
	cudaEventRecord(start_tDiv, 0);
#endif

  	
	// Launch the kernel
	kernel_trialDiv<<<dimGrid, dimBlock>>>(d_n1, d_r1);

	
#if PRINT_TIME
	// Stop and destroy the timer for MMM
	cudaEventRecord(stop_tDiv,0);
	cudaEventSynchronize(stop_tDiv);
	cudaEventElapsedTime(&elapsed_gpu, start_tDiv, stop_tDiv);
	printf("\nGPU time - kernel only: %f (msec)\n", elapsed_gpu);
	cudaEventDestroy(start_tDiv);
	cudaEventDestroy(stop_tDiv);
#endif
	
	
	// Check for errors during launch
	CUDA_SAFE_CALL(cudaPeekAtLastError());
	
	// Transfer the results back to the host
	CUDA_SAFE_CALL(cudaMemcpy(h_n1, d_n1, allocSizeLL, cudaMemcpyDeviceToHost));
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
	
	// Allocate memory for n2,r2 on host
	h_n2 = (long long int *) malloc(allocSizeLL);
	h_r2 = (int *) malloc(allocSizeInt);
	
	arrLen = pruneArray(h_n1,h_n2,h_r1,arrLen);						//copy h_n1 to h_n2, only potential primes. arrLen gets length of new array
	
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_n2, allocSizeLL));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_r2, allocSizeInt));	

	// Free-up host memory
	free(h_n1);
	free(h_r1);

// done with trial division	
////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
	
//	printf("potential primes: \n");	
//	for (i = 0; i < arrLen; i++)		
//		printf("%lld\t",h_n2[i]);


		
	free(h_n2);
	free(h_r2);	
	CUDA_SAFE_CALL(cudaFree(d_n2));
	CUDA_SAFE_CALL(cudaFree(d_r2));	
		
		
	return 0;

	
}

void initializeArray1D(long long int *arr, int len, long long int start) 		//each element increments by 2
{
	int i;
	for (i = 0; i < len; i++)
	{
		arr[i] = start;
		start += 2;
	}
}

void initializeArrayConst(int *arr, int len)				//initialize result matrix to 2
{
	int i;
	for (i = 0; i < len; i++)
		arr[i] = 1;
}

int pruneArray(long long int *arrSrc, long long int *arrDest, int *result, int len)
{
	int i;
	int length = 0;					//length of new array
	for (i = 0; i < len; i++)
	{
		if (result[i] == 1)							//only keep potential primes
		{
			arrDest[length] = arrSrc[i];			//write to the new array with current value
			length = length + 1;
		}	
	}
	return (length);
}