#ifndef CUDA_BULLETPROOF_H
#define CUDA_BULLETPROOF_H

#include "curve25519_ops.h"
#include "bulletproof_vectors.h"
#include "bulletproof_range_proof.h"

#ifdef __cplusplus
extern "C" {
#endif

// Multi-scalar multiplication optimization
void cuda_point_vector_multi_scalar_mul(ge25519* result,
                                      const FieldVector* scalars,
                                      const PointVector* points);

void cuda_point_vector_multi_scalar_mul_shared(ge25519* result,
                                             const FieldVector* scalars,
                                             const PointVector* points);

// Inner product calculation optimization
void cuda_field_vector_inner_product(fe25519* result,
                                   const FieldVector* a,
                                   const FieldVector* b);

void cuda_field_vector_inner_product_shared(fe25519* result,
                                          const FieldVector* a,
                                          const FieldVector* b);

// Field element operations optimization
void cuda_batch_field_add(fe25519* results,
                        const fe25519* a,
                        const fe25519* b,
                        size_t count);

void cuda_batch_field_sub(fe25519* results,
                        const fe25519* a,
                        const fe25519* b,
                        size_t count);

void cuda_batch_field_mul(fe25519* results,
                        const fe25519* a,
                        const fe25519* b,
                        size_t count);

void cuda_batch_field_square(fe25519* results,
                           const fe25519* inputs,
                           size_t count);

void cuda_batch_field_invert(fe25519* results,
                           const fe25519* inputs,
                           size_t count);

// Optimized SoA (Structure of Arrays) variants
void cuda_soa_field_add(fe25519* results,
                      const fe25519* a,
                      const fe25519* b,
                      size_t count);

// Optimized single proof verification
bool cuda_range_proof_verify(
    const RangeProof* proof,
    const ge25519* V,
    size_t n,
    const PointVector* G,
    const PointVector* H,
    const ge25519* g,
    const ge25519* h
);

// Optimized inner product verification
bool cuda_inner_product_verify(
    const InnerProductProof* proof,
    const ge25519* P,
    const PointVector* G,
    const PointVector* H,
    const ge25519* Q
);

// Performance benchmarking
void cuda_benchmark_multi_scalar_mul(int iterations, size_t vector_size);
void cuda_benchmark_inner_product(int iterations, size_t vector_size);
void cuda_benchmark_field_operations(int iterations, size_t batch_size);
void cuda_benchmark_range_proof(int iterations, size_t bit_size);

#ifdef __cplusplus
}
#endif

#endif // CUDA_BULLETPROOF_H
