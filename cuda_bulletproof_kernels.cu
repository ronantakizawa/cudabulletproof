
#include "curve25519_ops.h"
#include "bulletproof_vectors.h"
#include "device_curve25519_ops.cuh"  // Include the device operations
#include <stdio.h>

// Configuration parameters
#define BLOCK_SIZE 256
#define REDUCE_BLOCK_SIZE 128
#define MAX_SHARED_POINTS 64

// Common CUDA error checking macro
#define CUDA_CHECK(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, \
                    cudaGetErrorString(error)); \
            exit(EXIT_FAILURE); \
        } \
    } while(0)

/////////////// MULTI-SCALAR MULTIPLICATION OPTIMIZATION ///////////////

// Step 1: Kernel to perform individual scalar multiplications in parallel
__global__ void point_scalar_mul_kernel(ge25519* results,
                                       const fe25519* scalars,
                                       const ge25519* points,
                                       size_t n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        // Convert scalar to bytes for point multiplication
        uint8_t scalar_bytes[32];
        device_fe25519_tobytes(scalar_bytes, &scalars[idx]);

        // Perform scalar multiplication
        device_ge25519_scalarmult(&results[idx], scalar_bytes, &points[idx]);

        // Normalize result
        device_ge25519_normalize(&results[idx]);
    }
}

// Step 2: Kernel for parallel point addition with tree-structured reduction
__global__ void point_accumulate_kernel(ge25519* points,
                                      size_t n,
                                      size_t stride) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < n && idx + stride < n) {
        device_ge25519_add(&points[idx], &points[idx], &points[idx + stride]);
        device_ge25519_normalize(&points[idx]);
    }
}

// Forward declaration of helper function
void cudaPointVectorMultiScalarMulShared(ge25519* result,
                                       const FieldVector* scalars,
                                       const PointVector* points);

// Wrapper function to call the CUDA kernels
extern "C" void cuda_point_vector_multi_scalar_mul(ge25519* result,
                                                 const FieldVector* scalars,
                                                 const PointVector* points) {
    if (scalars->length != points->length) {
        fprintf(stderr, "Error: Vector lengths must match for multi-scalar multiplication\n");
        return;
    }

    size_t n = scalars->length;

    // Allocate device memory
    ge25519* d_points;
    fe25519* d_scalars;
    ge25519* d_results;

    CUDA_CHECK(cudaMalloc((void**)&d_points, n * sizeof(ge25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_scalars, n * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_results, n * sizeof(ge25519)));

    // Copy data to device
    CUDA_CHECK(cudaMemcpy(d_points, points->elements, n * sizeof(ge25519), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_scalars, scalars->elements, n * sizeof(fe25519), cudaMemcpyHostToDevice));

    // Calculate grid dimensions
    dim3 blockDim(BLOCK_SIZE);
    dim3 gridDim((n + BLOCK_SIZE - 1) / BLOCK_SIZE);

    // Step 1: Perform scalar multiplications in parallel
    point_scalar_mul_kernel<<<gridDim, blockDim>>>(d_results, d_scalars, d_points, n);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Step 2: Accumulate results using parallel reduction
    ge25519* d_temp;
    CUDA_CHECK(cudaMalloc((void**)&d_temp, n * sizeof(ge25519)));
    CUDA_CHECK(cudaMemcpy(d_temp, d_results, n * sizeof(ge25519), cudaMemcpyDeviceToDevice));

    for (size_t stride = 1; stride < n; stride *= 2) {
        size_t active_threads = n / (2 * stride);
        dim3 reduction_grid((active_threads + BLOCK_SIZE - 1) / BLOCK_SIZE);

        point_accumulate_kernel<<<reduction_grid, blockDim>>>(d_temp, n, stride);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());
    }

    // Copy the final result (the first element of d_temp)
    CUDA_CHECK(cudaMemcpy(result, d_temp, sizeof(ge25519), cudaMemcpyDeviceToHost));

    // Free device memory
    CUDA_CHECK(cudaFree(d_points));
    CUDA_CHECK(cudaFree(d_scalars));
    CUDA_CHECK(cudaFree(d_results));
    CUDA_CHECK(cudaFree(d_temp));
}

// Function using shared memory for better performance with small to medium sized inputs
extern "C" void cuda_point_vector_multi_scalar_mul_shared(ge25519* result,
                                                       const FieldVector* scalars,
                                                       const PointVector* points) {
    if (scalars->length != points->length) {
        fprintf(stderr, "Error: Vector lengths must match for multi-scalar multiplication\n");
        return;
    }

    size_t n = scalars->length;

    // For small n, use shared memory version
    if (n <= MAX_SHARED_POINTS) {
        // Implement optimized shared memory version
        cudaPointVectorMultiScalarMulShared(result, scalars, points);
        return;
    }

    // For larger n, use the standard version
    cuda_point_vector_multi_scalar_mul(result, scalars, points);
}

// Kernel using shared memory for small to medium inputs
__global__ void point_multi_scalar_mul_shared_kernel(ge25519* result,
                                                  const fe25519* scalars,
                                                  const ge25519* points,
                                                  size_t n) {
    __shared__ ge25519 shared_results[MAX_SHARED_POINTS];

    int tid = threadIdx.x;

    // Initialize shared memory
    if (tid < n) {
        // Convert scalar to bytes
        uint8_t scalar_bytes[32];
        device_fe25519_tobytes(scalar_bytes, &scalars[tid]);

        // Compute scalar multiplication
        device_ge25519_scalarmult(&shared_results[tid], scalar_bytes, &points[tid]);
        device_ge25519_normalize(&shared_results[tid]);
    }
    __syncthreads();

    // Parallel reduction in shared memory
    for (int stride = 1; stride < n; stride *= 2) {
        if (tid % (2 * stride) == 0 && tid + stride < n) {
            device_ge25519_add(&shared_results[tid], &shared_results[tid], &shared_results[tid + stride]);
            device_ge25519_normalize(&shared_results[tid]);
        }
        __syncthreads();
    }

    // Copy the final result
    if (tid == 0) {
        device_ge25519_copy(result, &shared_results[0]);
    }
}

// Wrapper for the shared memory kernel
void cudaPointVectorMultiScalarMulShared(ge25519* result,
                                       const FieldVector* scalars,
                                       const PointVector* points) {
    size_t n = scalars->length;

    // Allocate device memory
    ge25519* d_points;
    fe25519* d_scalars;
    ge25519* d_result;

    CUDA_CHECK(cudaMalloc((void**)&d_points, n * sizeof(ge25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_scalars, n * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_result, sizeof(ge25519)));

    // Copy data to device
    CUDA_CHECK(cudaMemcpy(d_points, points->elements, n * sizeof(ge25519), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_scalars, scalars->elements, n * sizeof(fe25519), cudaMemcpyHostToDevice));

    // Launch kernel with exact thread count needed
    point_multi_scalar_mul_shared_kernel<<<1, n>>>(d_result, d_scalars, d_points, n);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy result back
    CUDA_CHECK(cudaMemcpy(result, d_result, sizeof(ge25519), cudaMemcpyDeviceToHost));

    // Free device memory
    CUDA_CHECK(cudaFree(d_points));
    CUDA_CHECK(cudaFree(d_scalars));
    CUDA_CHECK(cudaFree(d_result));
}

// Warp-level implementation for even better performance - fixed for volatile issues
__device__ void warp_reduce_point(ge25519* result, ge25519* shared_points) {
    int lane = threadIdx.x & 31;  // Lane index within the warp

    // Perform reduction at warp level using warp shuffle
    for (int offset = 16; offset > 0; offset /= 2) {
        if (lane < offset) {
            // Make a non-volatile copy for operations
            ge25519 temp1 = shared_points[lane];
            ge25519 temp2 = shared_points[lane + offset];
            device_ge25519_add(&temp1, &temp1, &temp2);
            device_ge25519_normalize(&temp1);
            // Copy back to shared memory
            shared_points[lane] = temp1;
        }
        // Implicit warp synchronization (no __syncwarp needed for compute capability >= 7.0)
    }

    // First thread in the warp has the result
    if (lane == 0) {
        device_ge25519_copy(result, &shared_points[0]);
    }
}
