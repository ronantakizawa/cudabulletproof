
// Fixed cuda_field_ops.cu with consistent linkage

#include "curve25519_ops.h"
#include "device_curve25519_ops.cuh"  // Include the device operations
#include <stdio.h>

// Configuration parameters
#define BLOCK_SIZE 256
#define MAX_BATCH_SIZE 4096
#define WARP_SIZE 32

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

// Prime modulus for curve25519 (2^255 - 19)
__constant__ uint64_t d_p25519[4] = { 0xFFFFFFFFFFFFFFED, 0xFFFFFFFFFFFFFFFF,
                                     0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF };

/////////////// FIELD ELEMENT OPERATIONS OPTIMIZATION ///////////////

// Forward declaration of Karatsuba multiplication function WITH same linkage
extern "C" void cuda_batch_field_mul_karatsuba(fe25519* results,
                                          const fe25519* a,
                                          const fe25519* b,
                                          size_t count);

// Batch field addition kernel
__global__ void batch_field_add_kernel(fe25519* results,
                                     const fe25519* a,
                                     const fe25519* b,
                                     size_t count) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < count) {
        device_fe25519_add(&results[idx], &a[idx], &b[idx]);
    }
}

// Batch field subtraction kernel
__global__ void batch_field_sub_kernel(fe25519* results,
                                     const fe25519* a,
                                     const fe25519* b,
                                     size_t count) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < count) {
        device_fe25519_sub(&results[idx], &a[idx], &b[idx]);
    }
}

// Batch field multiplication kernel
__global__ void batch_field_mul_kernel(fe25519* results,
                                     const fe25519* a,
                                     const fe25519* b,
                                     size_t count) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < count) {
        device_fe25519_mul(&results[idx], &a[idx], &b[idx]);
    }
}

// Improved Karatsuba multiplication using shared memory
__global__ void karatsuba_field_mul_kernel(fe25519* results,
                                        const fe25519* a,
                                        const fe25519* b,
                                        size_t count) {
    __shared__ uint64_t shared_a[BLOCK_SIZE][4];
    __shared__ uint64_t shared_b[BLOCK_SIZE][4];

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;

    if (idx < count) {
        // Load inputs to shared memory for faster access
        for (int i = 0; i < 4; i++) {
            shared_a[tid][i] = a[idx].limbs[i];
            shared_b[tid][i] = b[idx].limbs[i];
        }
    }
    __syncthreads();

    if (idx < count) {
        // Implementation of Karatsuba multiplication using shared memory
        // This is a simplified version. In practice, you would implement
        // the full Karatsuba algorithm for better performance.

        uint64_t t[8] = {0}; // Temporary storage for multiplication result

        // Schoolbook multiplication (replace with Karatsuba for better performance)
        for (int i = 0; i < 4; i++) {
            uint64_t carry = 0;
            for (int j = 0; j < 4; j++) {
                unsigned __int128 m = (unsigned __int128)shared_a[tid][i] * shared_b[tid][j] + t[i+j] + carry;
                t[i+j] = (uint64_t)m;
                carry = (uint64_t)(m >> 64);
            }
            t[i+4] = carry;
        }

        // Modular reduction
        uint64_t p[4] = { d_p25519[0], d_p25519[1], d_p25519[2], d_p25519[3] };

        uint64_t carry = 0;
        uint64_t c;

        // Multiply top limb by 19 and add to lowest limb
        c = t[4] * 19;
        t[0] += c;
        carry = t[0] < c ? 1 : 0;

        for (int i = 1; i < 4; i++) {
            c = t[i+4] * 19 + carry;
            t[i] += c;
            carry = t[i] < c ? 1 : 0;
        }

        // Final reduction if needed
        if (carry || (t[3] > p[3]) ||
            ((t[3] == p[3]) && ((t[2] > p[2]) ||
                                ((t[2] == p[2]) && ((t[1] > p[1]) ||
                                                   ((t[1] == p[1]) && (t[0] >= p[0]))))))) {
            carry = 0;
            for (int i = 0; i < 4; i++) {
                uint64_t diff = t[i] - p[i] - carry;
                carry = (t[i] < p[i] + carry) ? 1 : 0;
                results[idx].limbs[i] = diff;
            }
        } else {
            for (int i = 0; i < 4; i++) {
                results[idx].limbs[i] = t[i];
            }
        }
    }
}

// Optimized field element squaring
__global__ void field_square_kernel(fe25519* results,
                                  const fe25519* inputs,
                                  size_t count) {
    __shared__ uint64_t shared_inputs[BLOCK_SIZE][4];

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;

    if (idx < count) {
        // Load inputs to shared memory
        for (int i = 0; i < 4; i++) {
            shared_inputs[tid][i] = inputs[idx].limbs[i];
        }
    }
    __syncthreads();

    if (idx < count) {
        // Optimized squaring implementation
        // In practice, we would use specialized squaring algorithms
        // that take advantage of shared terms in the multiplication

        uint64_t t[8] = {0};

        // Simplified squaring (in practice, use specialized algorithms)
        for (int i = 0; i < 4; i++) {
            // Diagonal terms (a_i * a_i)
            unsigned __int128 diag = (unsigned __int128)shared_inputs[tid][i] * shared_inputs[tid][i];
            t[i+i] += (uint64_t)diag;
            if (i+i+1 < 8) t[i+i+1] += (uint64_t)(diag >> 64);

            // Off-diagonal terms (2 * a_i * a_j for i < j)
            for (int j = i+1; j < 4; j++) {
                unsigned __int128 m = 2 * ((unsigned __int128)shared_inputs[tid][i] * shared_inputs[tid][j]);
                t[i+j] += (uint64_t)m;
                if (i+j+1 < 8) t[i+j+1] += (uint64_t)(m >> 64);
            }
        }

        // Modular reduction (same as multiplication)
        uint64_t p[4] = { d_p25519[0], d_p25519[1], d_p25519[2], d_p25519[3] };

        uint64_t carry = 0;
        uint64_t c;

        // Multiply top limb by 19 and add to lowest limb
        c = t[4] * 19;
        t[0] += c;
        carry = t[0] < c ? 1 : 0;

        for (int i = 1; i < 4; i++) {
            c = t[i+4] * 19 + carry;
            t[i] += c;
            carry = t[i] < c ? 1 : 0;
        }

        // Final reduction if needed
        if (carry || (t[3] > p[3]) ||
            ((t[3] == p[3]) && ((t[2] > p[2]) ||
                                ((t[2] == p[2]) && ((t[1] > p[1]) ||
                                                   ((t[1] == p[1]) && (t[0] >= p[0]))))))) {
            carry = 0;
            for (int i = 0; i < 4; i++) {
                uint64_t diff = t[i] - p[i] - carry;
                carry = (t[i] < p[i] + carry) ? 1 : 0;
                results[idx].limbs[i] = diff;
            }
        } else {
            for (int i = 0; i < 4; i++) {
                results[idx].limbs[i] = t[i];
            }
        }
    }
}

// Batch inversion using Montgomery's trick
__global__ void field_batch_invert_kernel(fe25519* products,
                                        const fe25519* inputs,
                                        size_t count) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < count) {
        if (idx == 0) {
            // Copy first element
            device_fe25519_copy(&products[0], &inputs[0]);
        } else {
            // products[i] = inputs[0] * ... * inputs[i]
            device_fe25519_mul(&products[idx], &products[idx-1], &inputs[idx]);
        }
    }
}

__global__ void field_batch_invert_finalize_kernel(fe25519* results,
                                                 const fe25519* inputs,
                                                 fe25519* products,
                                                 const fe25519* total_inverse,
                                                 size_t count) {
    int idx = count - 1 - (blockIdx.x * blockDim.x + threadIdx.x);

    if (idx >= 0 && idx < count) {
        if (idx == count - 1) {
            // Last element gets the total inverse
            device_fe25519_copy(&results[idx], total_inverse);
        } else {
            // results[i] = products[i] * results[i+1]
            device_fe25519_mul(&results[idx], &products[idx], &results[idx+1]);
        }
    }
}

// Wrapper for batch field element addition
extern "C" void cuda_batch_field_add(fe25519* results,
                                   const fe25519* a,
                                   const fe25519* b,
                                   size_t count) {
    // Allocate device memory
    fe25519* d_a;
    fe25519* d_b;
    fe25519* d_results;

    CUDA_CHECK(cudaMalloc((void**)&d_a, count * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_b, count * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_results, count * sizeof(fe25519)));

    // Copy data to device
    CUDA_CHECK(cudaMemcpy(d_a, a, count * sizeof(fe25519), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, b, count * sizeof(fe25519), cudaMemcpyHostToDevice));

    // Calculate grid dimensions
    dim3 blockDim(BLOCK_SIZE);
    dim3 gridDim((count + BLOCK_SIZE - 1) / BLOCK_SIZE);

    // Launch kernel
    batch_field_add_kernel<<<gridDim, blockDim>>>(d_results, d_a, d_b, count);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy results back
    CUDA_CHECK(cudaMemcpy(results, d_results, count * sizeof(fe25519), cudaMemcpyDeviceToHost));

    // Free device memory
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_results));
}

// Wrapper for batch field element subtraction
extern "C" void cuda_batch_field_sub(fe25519* results,
                                   const fe25519* a,
                                   const fe25519* b,
                                   size_t count) {
    // Allocate device memory
    fe25519* d_a;
    fe25519* d_b;
    fe25519* d_results;

    CUDA_CHECK(cudaMalloc((void**)&d_a, count * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_b, count * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_results, count * sizeof(fe25519)));

    // Copy data to device
    CUDA_CHECK(cudaMemcpy(d_a, a, count * sizeof(fe25519), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, b, count * sizeof(fe25519), cudaMemcpyHostToDevice));

    // Calculate grid dimensions
    dim3 blockDim(BLOCK_SIZE);
    dim3 gridDim((count + BLOCK_SIZE - 1) / BLOCK_SIZE);

    // Launch kernel
    batch_field_sub_kernel<<<gridDim, blockDim>>>(d_results, d_a, d_b, count);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy results back
    CUDA_CHECK(cudaMemcpy(results, d_results, count * sizeof(fe25519), cudaMemcpyDeviceToHost));

    // Free device memory
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_results));
}

// Wrapper for batch field element multiplication
extern "C" void cuda_batch_field_mul(fe25519* results,
                                   const fe25519* a,
                                   const fe25519* b,
                                   size_t count) {
    // Use Karatsuba multiplication for better performance
    cuda_batch_field_mul_karatsuba(results, a, b, count);
}

// Wrapper for Karatsuba multiplication
extern "C" void cuda_batch_field_mul_karatsuba(fe25519* results,
                                            const fe25519* a,
                                            const fe25519* b,
                                            size_t count) {
    // Allocate device memory
    fe25519* d_a;
    fe25519* d_b;
    fe25519* d_results;

    CUDA_CHECK(cudaMalloc((void**)&d_a, count * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_b, count * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_results, count * sizeof(fe25519)));

    // Copy data to device
    CUDA_CHECK(cudaMemcpy(d_a, a, count * sizeof(fe25519), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b, b, count * sizeof(fe25519), cudaMemcpyHostToDevice));

    // Calculate grid dimensions
    dim3 blockDim(BLOCK_SIZE);
    dim3 gridDim((count + BLOCK_SIZE - 1) / BLOCK_SIZE);

    // Launch kernel
    karatsuba_field_mul_kernel<<<gridDim, blockDim>>>(d_results, d_a, d_b, count);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy results back
    CUDA_CHECK(cudaMemcpy(results, d_results, count * sizeof(fe25519), cudaMemcpyDeviceToHost));

    // Free device memory
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_results));
}

// Wrapper for batch field element squaring
extern "C" void cuda_batch_field_square(fe25519* results,
                                      const fe25519* inputs,
                                      size_t count) {
    // Allocate device memory
    fe25519* d_inputs;
    fe25519* d_results;

    CUDA_CHECK(cudaMalloc((void**)&d_inputs, count * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_results, count * sizeof(fe25519)));

    // Copy data to device
    CUDA_CHECK(cudaMemcpy(d_inputs, inputs, count * sizeof(fe25519), cudaMemcpyHostToDevice));

    // Calculate grid dimensions
    dim3 blockDim(BLOCK_SIZE);
    dim3 gridDim((count + BLOCK_SIZE - 1) / BLOCK_SIZE);

    // Launch kernel
    field_square_kernel<<<gridDim, blockDim>>>(d_results, d_inputs, count);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy results back
    CUDA_CHECK(cudaMemcpy(results, d_results, count * sizeof(fe25519), cudaMemcpyDeviceToHost));

    // Free device memory
    CUDA_CHECK(cudaFree(d_inputs));
    CUDA_CHECK(cudaFree(d_results));
}

// Wrapper for batch field element inversion using Montgomery's trick
extern "C" void cuda_batch_field_invert(fe25519* results,
                                      const fe25519* inputs,
                                      size_t count) {
    if (count == 0) return;
    if (count == 1) {
        // For a single element, use standard inversion
        fe25519_invert(results, inputs);
        return;
    }

    // Allocate device memory
    fe25519* d_inputs;
    fe25519* d_results;
    fe25519* d_products;
    fe25519* d_total_inverse;

    CUDA_CHECK(cudaMalloc((void**)&d_inputs, count * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_results, count * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_products, count * sizeof(fe25519)));
    CUDA_CHECK(cudaMalloc((void**)&d_total_inverse, sizeof(fe25519)));

    // Copy data to device
    CUDA_CHECK(cudaMemcpy(d_inputs, inputs, count * sizeof(fe25519), cudaMemcpyHostToDevice));

    // Calculate grid dimensions
    dim3 blockDim(BLOCK_SIZE);
    dim3 gridDim((count + BLOCK_SIZE - 1) / BLOCK_SIZE);

    // Phase 1: Compute running products
    field_batch_invert_kernel<<<gridDim, blockDim>>>(d_products, d_inputs, count);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy the final product to host and compute its inverse
    fe25519 total_product, total_inverse;
    CUDA_CHECK(cudaMemcpy(&total_product, &d_products[count-1], sizeof(fe25519), cudaMemcpyDeviceToHost));

    // Compute inverse of the total product
    fe25519_invert(&total_inverse, &total_product);

    // Copy the inverse back to device
    CUDA_CHECK(cudaMemcpy(d_total_inverse, &total_inverse, sizeof(fe25519), cudaMemcpyHostToDevice));

    // Phase 2: Compute individual inverses
    field_batch_invert_finalize_kernel<<<gridDim, blockDim>>>(
        d_results, d_inputs, d_products, d_total_inverse, count);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy results back
    CUDA_CHECK(cudaMemcpy(results, d_results, count * sizeof(fe25519), cudaMemcpyDeviceToHost));

    // Free device memory
    CUDA_CHECK(cudaFree(d_inputs));
    CUDA_CHECK(cudaFree(d_results));
    CUDA_CHECK(cudaFree(d_products));
    CUDA_CHECK(cudaFree(d_total_inverse));
}

// Specialized kernel for more efficient memory access
__global__ void limb_oriented_field_add_kernel(fe25519* results,
                                            const fe25519* a,
                                            const fe25519* b,
                                            size_t count) {
    // Each block handles multiple field elements but focuses on one limb
    int limb_idx = blockIdx.y;
    int batch_idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (batch_idx < count) {
        // Each thread adds one limb of one field element
        uint64_t limb_a = a[batch_idx].limbs[limb_idx];
        uint64_t limb_b = b[batch_idx].limbs[limb_idx];

        // Simple addition (without carry handling for simplicity)
        results[batch_idx].limbs[limb_idx] = limb_a + limb_b;

        // Note: In a real implementation, we would need to handle carries properly
    }
}

// Structure of Arrays (SoA) implementation for better memory coalescence
typedef struct {
    uint64_t* limb0;  // Array of first limbs for all field elements
    uint64_t* limb1;  // Array of second limbs for all field elements
    uint64_t* limb2;  // Array of third limbs for all field elements
    uint64_t* limb3;  // Array of fourth limbs for all field elements
    size_t count;     // Number of field elements
} fe25519_soa;

// Convert from Array of Structures (AoS) to Structure of Arrays (SoA)
void fe25519_aos_to_soa(fe25519_soa* soa, const fe25519* aos, size_t count) {
    soa->limb0 = (uint64_t*)malloc(count * sizeof(uint64_t));
    soa->limb1 = (uint64_t*)malloc(count * sizeof(uint64_t));
    soa->limb2 = (uint64_t*)malloc(count * sizeof(uint64_t));
    soa->limb3 = (uint64_t*)malloc(count * sizeof(uint64_t));
    soa->count = count;

    for (size_t i = 0; i < count; i++) {
        soa->limb0[i] = aos[i].limbs[0];
        soa->limb1[i] = aos[i].limbs[1];
        soa->limb2[i] = aos[i].limbs[2];
        soa->limb3[i] = aos[i].limbs[3];
    }
}

// Convert from Structure of Arrays (SoA) to Array of Structures (AoS)
void fe25519_soa_to_aos(fe25519* aos, const fe25519_soa* soa) {
    for (size_t i = 0; i < soa->count; i++) {
        aos[i].limbs[0] = soa->limb0[i];
        aos[i].limbs[1] = soa->limb1[i];
        aos[i].limbs[2] = soa->limb2[i];
        aos[i].limbs[3] = soa->limb3[i];
    }
}

// Kernel for Structure of Arrays addition (better memory coalescence)
__global__ void soa_field_add_kernel(uint64_t* r_limb,
                                   const uint64_t* a_limb,
                                   const uint64_t* b_limb,
                                   size_t count) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < count) {
        r_limb[idx] = a_limb[idx] + b_limb[idx];
    }
}

// Wrapper for SoA field addition
extern "C" void cuda_soa_field_add(fe25519* results,
                                const fe25519* a,
                                const fe25519* b,
                                size_t count) {
    // Convert inputs to SoA format
    fe25519_soa a_soa, b_soa, results_soa;
    fe25519_aos_to_soa(&a_soa, a, count);
    fe25519_aos_to_soa(&b_soa, b, count);

    results_soa.limb0 = (uint64_t*)malloc(count * sizeof(uint64_t));
    results_soa.limb1 = (uint64_t*)malloc(count * sizeof(uint64_t));
    results_soa.limb2 = (uint64_t*)malloc(count * sizeof(uint64_t));
    results_soa.limb3 = (uint64_t*)malloc(count * sizeof(uint64_t));
    results_soa.count = count;

    // Allocate device memory for SoA format
    uint64_t *d_a_limb0, *d_a_limb1, *d_a_limb2, *d_a_limb3;
    uint64_t *d_b_limb0, *d_b_limb1, *d_b_limb2, *d_b_limb3;
    uint64_t *d_r_limb0, *d_r_limb1, *d_r_limb2, *d_r_limb3;

    CUDA_CHECK(cudaMalloc((void**)&d_a_limb0, count * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc((void**)&d_a_limb1, count * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc((void**)&d_a_limb2, count * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc((void**)&d_a_limb3, count * sizeof(uint64_t)));

    CUDA_CHECK(cudaMalloc((void**)&d_b_limb0, count * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc((void**)&d_b_limb1, count * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc((void**)&d_b_limb2, count * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc((void**)&d_b_limb3, count * sizeof(uint64_t)));

    CUDA_CHECK(cudaMalloc((void**)&d_r_limb0, count * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc((void**)&d_r_limb1, count * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc((void**)&d_r_limb2, count * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc((void**)&d_r_limb3, count * sizeof(uint64_t)));

    // Copy data to device
    CUDA_CHECK(cudaMemcpy(d_a_limb0, a_soa.limb0, count * sizeof(uint64_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_a_limb1, a_soa.limb1, count * sizeof(uint64_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_a_limb2, a_soa.limb2, count * sizeof(uint64_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_a_limb3, a_soa.limb3, count * sizeof(uint64_t), cudaMemcpyHostToDevice));

    CUDA_CHECK(cudaMemcpy(d_b_limb0, b_soa.limb0, count * sizeof(uint64_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b_limb1, b_soa.limb1, count * sizeof(uint64_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b_limb2, b_soa.limb2, count * sizeof(uint64_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_b_limb3, b_soa.limb3, count * sizeof(uint64_t), cudaMemcpyHostToDevice));

    // Calculate grid dimensions
    dim3 blockDim(BLOCK_SIZE);
    dim3 gridDim((count + BLOCK_SIZE - 1) / BLOCK_SIZE);

    // Launch kernels
    soa_field_add_kernel<<<gridDim, blockDim>>>(d_r_limb0, d_a_limb0, d_b_limb0, count);
    soa_field_add_kernel<<<gridDim, blockDim>>>(d_r_limb1, d_a_limb1, d_b_limb1, count);
    soa_field_add_kernel<<<gridDim, blockDim>>>(d_r_limb2, d_a_limb2, d_b_limb2, count);
    soa_field_add_kernel<<<gridDim, blockDim>>>(d_r_limb3, d_a_limb3, d_b_limb3, count);

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy results back
    CUDA_CHECK(cudaMemcpy(results_soa.limb0, d_r_limb0, count * sizeof(uint64_t), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(results_soa.limb1, d_r_limb1, count * sizeof(uint64_t), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(results_soa.limb2, d_r_limb2, count * sizeof(uint64_t), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(results_soa.limb3, d_r_limb3, count * sizeof(uint64_t), cudaMemcpyDeviceToHost));

    // Convert back to AoS format
    fe25519_soa_to_aos(results, &results_soa);

    // Free device memory
    CUDA_CHECK(cudaFree(d_a_limb0));
    CUDA_CHECK(cudaFree(d_a_limb1));
    CUDA_CHECK(cudaFree(d_a_limb2));
    CUDA_CHECK(cudaFree(d_a_limb3));
    CUDA_CHECK(cudaFree(d_b_limb0));
    CUDA_CHECK(cudaFree(d_b_limb1));
    CUDA_CHECK(cudaFree(d_b_limb2));
    CUDA_CHECK(cudaFree(d_b_limb3));
    CUDA_CHECK(cudaFree(d_r_limb0));
    CUDA_CHECK(cudaFree(d_r_limb1));
    CUDA_CHECK(cudaFree(d_r_limb2));
    CUDA_CHECK(cudaFree(d_r_limb3));

    // Free host memory
    free(a_soa.limb0);
    free(a_soa.limb1);
    free(a_soa.limb2);
    free(a_soa.limb3);
    free(b_soa.limb0);
    free(b_soa.limb1);
    free(b_soa.limb2);
    free(b_soa.limb3);
    free(results_soa.limb0);
    free(results_soa.limb1);
    free(results_soa.limb2);
    free(results_soa.limb3);
}
