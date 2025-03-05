
// Fixed cuda_inner_product.cu

#include "curve25519_ops.h"
#include "bulletproof_vectors.h"
#include "device_curve25519_ops.cuh"  // Include the device operations
#include <stdio.h>

// Configuration parameters
#define BLOCK_SIZE 256
#define WARP_SIZE 32
#define MAX_SHARED_ELEMENTS 512

// CUDA error checking macro (reusing from previous file)
#define CUDA_CHECK(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, \
                    cudaGetErrorString(error)); \
            exit(EXIT_FAILURE); \
        } \
    } while(0)

/////////////// INNER PRODUCT CALCULATION OPTIMIZATION ///////////////

// Forward declaration of shared memory helper
void cuda_field_vector_inner_product_shared(fe25519* result,
                                          const FieldVector* a,
                                          const FieldVector* b);

// Kernel for computing inner products using parallel reduction
__global__ void field_vector_inner_product_kernel(fe25519* result,
                                                const fe25519* a,
                                                const fe25519* b,
                                                size_t n) {
    __shared__ fe25519 partial_sums[BLOCK_SIZE];

    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Initialize shared memory
    device_fe25519_0(&partial_sums[tid]);

    // Each thread computes one or more products
    while (idx < n) {
        fe25519 product;
        device_fe25519_mul(&product, &a[idx], &b[idx]);
        device_fe25519_add(&partial_sums[tid], &partial_sums[tid], &product);
        idx += gridDim.x * blockDim.x;
    }
    __syncthreads();

    // Perform reduction in shared memory
    for (int stride = BLOCK_SIZE/2; stride > 0; stride >>= 1) {
        if (tid < stride) {
            device_fe25519_add(&partial_sums[tid], &partial_sums[tid], &partial_sums[tid + stride]);
        }
        __syncthreads();
    }

    // Thread 0 writes the final result
    if (tid == 0) {
        device_fe25519_copy(result + blockIdx.x, &partial_sums[0]);
    }
}

// Second-level reduction kernel
__global__ void fe25519_reduce_kernel(fe25519* result, const fe25519* partial_results, size_t n) {
    __shared__ fe25519 shared_data[BLOCK_SIZE];

    int tid = threadIdx.x;

    // Copy to shared memory
    if (tid < n) {
        device_fe25519_copy(&shared_data[tid], &partial_results[tid]);
    } else {
        device_fe25519_0(&shared_data[tid]);
    }
    __syncthreads();

    // Reduction in shared memory
    for (int stride = BLOCK_SIZE/2; stride > 0; stride >>= 1) {
        if (tid < stride && tid + stride < n) {
            device_fe25519_add(&shared_data[tid], &shared_data[tid], &shared_data[tid + stride]);
        }
        __syncthreads();
    }

    // Thread 0 writes the final result
    if (tid == 0) {
        device_fe25519_copy(result, &shared_data[0]);
    }
}

// Optimized wrapper for field vector inner product
extern "C" void cuda_field_vector_inner_product(fe25519* result,
                                             const FieldVector* a,
                                             const FieldVector* b) {
    if (a->length != b->length) {
        fprintf(stderr, "Error: Vector lengths must match for inner product\n");
        return;
    }

    size_t n = a->length;

    // Use optimized shared memory version for small inputs
    if (n <= MAX_SHARED_ELEMENTS) {
        cuda_field_vector_inner_product_shared(result, a, b);
        return;
    }

    // Allocate device memory
    fe25519* d_a;
    fe25519* d_b;
    fe25519* d_result;
    fe25519* d_temp;

    CUDA_CHECK(cudaMalloc((void**)&d_a, n * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_b, n * sizeof(fe25519)));

    // Calculate grid dimensions
    int num_blocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
    num_blocks = (num_blocks > 1024) ? 1024 : num_blocks; // Limit max blocks

    CUDA_CHECK(cudaMalloc((void**)&d_temp, num_blocks * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_result, sizeof(fe25519)));

    // Copy data to device
    CUDA_CHECK(cudaMemcpy(d_a, a->elements, n * sizeof(fe25519), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, b->elements, n * sizeof(fe25519), cudaMemcpyHostToDevice));

    // Step 1: Compute partial inner products
    field_vector_inner_product_kernel<<<num_blocks, BLOCK_SIZE>>>(d_temp, d_a, d_b, n);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Step 2: Reduce partial results
    fe25519_reduce_kernel<<<1, BLOCK_SIZE>>>(d_result, d_temp, num_blocks);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy result back to host
    CUDA_CHECK(cudaMemcpy(result, d_result, sizeof(fe25519), cudaMemcpyDeviceToHost));

    // Free device memory
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_temp));
    CUDA_CHECK(cudaFree(d_result));
}

// Kernel for small inputs using shared memory
__global__ void field_vector_inner_product_shared_kernel(fe25519* result,
                                                      const fe25519* a,
                                                      const fe25519* b,
                                                      size_t n) {
    __shared__ fe25519 partial_products[MAX_SHARED_ELEMENTS];

    int tid = threadIdx.x;

    // Compute individual products
    if (tid < n) {
        device_fe25519_mul(&partial_products[tid], &a[tid], &b[tid]);
    } else {
        device_fe25519_0(&partial_products[tid]);
    }
    __syncthreads();

    // Parallel reduction within thread block
    for (unsigned int stride = blockDim.x/2; stride > 0; stride >>= 1) {
        if (tid < stride && tid + stride < n) {
            device_fe25519_add(&partial_products[tid], &partial_products[tid], &partial_products[tid + stride]);
        }
        __syncthreads();
    }

    // First thread writes the result
    if (tid == 0) {
        device_fe25519_copy(result, &partial_products[0]);
    }
}

// Wrapper for shared memory version
void cuda_field_vector_inner_product_shared(fe25519* result,
                                         const FieldVector* a,
                                         const FieldVector* b) {
    size_t n = a->length;

    // Allocate device memory
    fe25519* d_a;
    fe25519* d_b;
    fe25519* d_result;

    CUDA_CHECK(cudaMalloc((void**)&d_a, n * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_b, n * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_result, sizeof(fe25519)));

    // Copy data to device
    CUDA_CHECK(cudaMemcpy(d_a, a->elements, n * sizeof(fe25519), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, b->elements, n * sizeof(fe25519), cudaMemcpyHostToDevice));

    // Launch kernel with sufficient threads to cover the input size
    int num_threads = min((int)n, MAX_SHARED_ELEMENTS);
    field_vector_inner_product_shared_kernel<<<1, num_threads>>>(d_result, d_a, d_b, n);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy result back
    CUDA_CHECK(cudaMemcpy(result, d_result, sizeof(fe25519), cudaMemcpyDeviceToHost));

    // Free device memory
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_result));
}

// Warp-level primitives for faster reduction - fixed for volatile issues
__device__ void warp_reduce_field_element(fe25519* out, fe25519* sdata) {
    int lane = threadIdx.x & (WARP_SIZE - 1);

    // Use non-volatile copies to perform operations
    if (lane < 16) {
        fe25519 temp1 = sdata[lane];
        fe25519 temp2 = sdata[lane + 16];
        device_fe25519_add(&temp1, &temp1, &temp2);
        sdata[lane] = temp1;
    }
    if (lane < 8) {
        fe25519 temp1 = sdata[lane];
        fe25519 temp2 = sdata[lane + 8];
        device_fe25519_add(&temp1, &temp1, &temp2);
        sdata[lane] = temp1;
    }
    if (lane < 4) {
        fe25519 temp1 = sdata[lane];
        fe25519 temp2 = sdata[lane + 4];
        device_fe25519_add(&temp1, &temp1, &temp2);
        sdata[lane] = temp1;
    }
    if (lane < 2) {
        fe25519 temp1 = sdata[lane];
        fe25519 temp2 = sdata[lane + 2];
        device_fe25519_add(&temp1, &temp1, &temp2);
        sdata[lane] = temp1;
    }
    if (lane < 1) {
        fe25519 temp1 = sdata[lane];
        fe25519 temp2 = sdata[lane + 1];
        device_fe25519_add(&temp1, &temp1, &temp2);
        sdata[lane] = temp1;
    }

    if (lane == 0) {
        device_fe25519_copy(out, &sdata[0]);
    }
}

// Batch processing for multiple inner products
__global__ void batch_inner_product_kernel(fe25519* results,
                                        const fe25519* a_vectors,
                                        const fe25519* b_vectors,
                                        size_t n,
                                        size_t batch_size) {
    __shared__ fe25519 shared_sums[BLOCK_SIZE];

    int batch_id = blockIdx.y;
    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Base offsets for this batch
    const fe25519* a = a_vectors + batch_id * n;
    const fe25519* b = b_vectors + batch_id * n;

    // Initialize shared memory
    device_fe25519_0(&shared_sums[tid]);

    // Each thread computes one or more products
    while (idx < n) {
        fe25519 product;
        device_fe25519_mul(&product, &a[idx], &b[idx]);
        device_fe25519_add(&shared_sums[tid], &shared_sums[tid], &product);
        idx += gridDim.x * blockDim.x;
    }
    __syncthreads();

    // Perform reduction in shared memory
    for (int stride = BLOCK_SIZE/2; stride >= WARP_SIZE; stride >>= 1) {
        if (tid < stride) {
            device_fe25519_add(&shared_sums[tid], &shared_sums[tid], &shared_sums[tid + stride]);
        }
        __syncthreads();
    }

    // Final warp reduction
    if (tid < WARP_SIZE) {
        warp_reduce_field_element(&results[batch_id], shared_sums);
    }
}

// Wrapper for batch inner product calculation
extern "C" void cuda_batch_field_vector_inner_product(fe25519* results,
                                                   const FieldVector* a_vectors,
                                                   const FieldVector* b_vectors,
                                                   size_t num_vectors) {
    size_t n = a_vectors[0].length;

    // Prepare data in contiguous arrays
    fe25519* h_a_contiguous = (fe25519*)malloc(num_vectors * n * sizeof(fe25519));
    fe25519* h_b_contiguous = (fe25519*)malloc(num_vectors * n * sizeof(fe25519));

    for (size_t i = 0; i < num_vectors; i++) {
        memcpy(h_a_contiguous + i * n, a_vectors[i].elements, n * sizeof(fe25519));
        memcpy(h_b_contiguous + i * n, b_vectors[i].elements, n * sizeof(fe25519));
    }

    // Allocate device memory
    fe25519* d_a;
    fe25519* d_b;
    fe25519* d_results;

    CUDA_CHECK(cudaMalloc((void**)&d_a, num_vectors * n * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_b, num_vectors * n * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_results, num_vectors * sizeof(fe25519)));

    // Copy data to device
    CUDA_CHECK(cudaMemcpy(d_a, h_a_contiguous, num_vectors * n * sizeof(fe25519), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, h_b_contiguous, num_vectors * n * sizeof(fe25519), cudaMemcpyHostToDevice));

    // Configure grid
    dim3 block(BLOCK_SIZE);
    dim3 grid(min(1024, (int)((n + BLOCK_SIZE - 1) / BLOCK_SIZE)), num_vectors);

    // Launch kernel
    batch_inner_product_kernel<<<grid, block>>>(d_results, d_a, d_b, n, num_vectors);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy results back
    CUDA_CHECK(cudaMemcpy(results, d_results, num_vectors * sizeof(fe25519), cudaMemcpyDeviceToHost));

    // Free memory
    free(h_a_contiguous);
    free(h_b_contiguous);
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_results));
}
