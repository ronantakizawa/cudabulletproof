#ifndef BULLETPROOF_VECTORS_H
#define BULLETPROOF_VECTORS_H

#include "curve25519_ops.h"
#include <stdint.h>

// Structure to represent a vector of field elements
typedef struct {
    fe25519* elements;
    size_t length;
} FieldVector;

// Structure to represent a vector of points
typedef struct {
    ge25519* elements;
    size_t length;
} PointVector;

// Initialize a field vector of given size
void field_vector_init(FieldVector* vec, size_t length);

// Free memory allocated for field vector
void field_vector_free(FieldVector* vec);

// Set all elements to 0
void field_vector_clear(FieldVector* vec);

// Copy vector: dest = src
void field_vector_copy(FieldVector* dest, const FieldVector* src);

// Vector-scalar multiplication: result = scalar * vec
void field_vector_scalar_mul(FieldVector* result, const FieldVector* vec, const fe25519* scalar);

// Vector addition: result = a + b
void field_vector_add(FieldVector* result, const FieldVector* a, const FieldVector* b);

// Vector subtraction: result = a - b
void field_vector_sub(FieldVector* result, const FieldVector* a, const FieldVector* b);

// Vector inner product: result = <a, b>
void field_vector_inner_product(fe25519* result, const FieldVector* a, const FieldVector* b);

// Hadamard product: result = a ○ b (element-wise multiplication)
void field_vector_hadamard(FieldVector* result, const FieldVector* a, const FieldVector* b);

// Initialize a point vector of given size
void point_vector_init(PointVector* vec, size_t length);

// Free memory allocated for point vector
void point_vector_free(PointVector* vec);

// Set all elements to identity
void point_vector_clear(PointVector* vec);

// Copy vector: dest = src
void point_vector_copy(PointVector* dest, const PointVector* src);

// Vector-scalar multiplication: result = scalar * vec
void point_vector_scalar_mul(PointVector* result, const PointVector* vec, const fe25519* scalar);

// Multi-scalar multiplication: result = <scalars, points>
void point_vector_multi_scalar_mul(ge25519* result, const FieldVector* scalars, const PointVector* points);

// Inner product protocol structure
typedef struct {
    size_t n;               // Size of original vectors (power of 2)
    FieldVector a;          // Left vector
    FieldVector b;          // Right vector
    fe25519 c;              // Inner product value <a,b>
    PointVector L;          // Left commitments
    PointVector R;          // Right commitments
    size_t L_len;           // Length of L and R (log n)
    fe25519 x;              // Challenge
} InnerProductProof;

// Initialize an inner product proof
void inner_product_proof_init(InnerProductProof* proof, size_t n);

// Free memory allocated for an inner product proof
void inner_product_proof_free(InnerProductProof* proof);

/**
 * Generate an inner product proof
 *
 * @param proof Output parameter to store the generated proof
 * @param a_in Left vector for inner product
 * @param b_in Right vector for inner product
 * @param G Base point vector for left vector commitment
 * @param H Base point vector for right vector commitment
 * @param Q Additional base point for cross-term
 * @param c_in Claimed inner product value
 * @param transcript_hash Initial transcript state for Fiat-Shamir
 */
void inner_product_prove(
    InnerProductProof* proof,
    const FieldVector* a_in,
    const FieldVector* b_in,
    const PointVector* G,
    const PointVector* H,
    const ge25519* Q,
    const fe25519* c_in,
    const uint8_t* transcript_hash
);

/**
 * Verify an inner product proof
 *
 * @param proof The inner product proof to verify
 * @param P Point representing the commitment to verify against
 * @param G Base point vector for left vector commitment
 * @param H Base point vector for right vector commitment
 * @param Q Additional base point for cross-term
 * @return true if the proof is valid, false otherwise
 */
bool inner_product_verify(
    const InnerProductProof* proof,
    const ge25519* P,
    const PointVector* G,
    const PointVector* H,
    const ge25519* Q
);

#endif // BULLETPROOF_VECTORS_H
