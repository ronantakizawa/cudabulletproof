
#ifndef BULLETPROOF_RANGE_PROOF_H
#define BULLETPROOF_RANGE_PROOF_H

#include "curve25519_ops.h"
#include "bulletproof_vectors.h"

// Structure for a Bulletproof range proof
typedef struct {
    ge25519 V;           // Value commitment
    ge25519 A;           // Polynomial commitment for a
    ge25519 S;           // Polynomial commitment for s
    ge25519 T1;          // Polynomial commitment t1
    ge25519 T2;          // Polynomial commitment t2
    fe25519 taux;        // Blinding factor for t
    fe25519 mu;          // Blinding factor for inner product
    fe25519 t;           // Polynomial evaluation
    InnerProductProof ip_proof;  // Inner product proof
} RangeProof;

// Initialize a range proof
void range_proof_init(RangeProof* proof, size_t n);

// Free memory allocated for a range proof
void range_proof_free(RangeProof* proof);

// Helper function to generate a Pedersen commitment
void pedersen_commit(ge25519* result, const fe25519* value, const fe25519* blinding, const ge25519* g, const ge25519* h);

// Helper function to generate a vector of powers of a base value
void powers_of(FieldVector* result, const fe25519* base, size_t n);

// Compute precise delta value for polynomial identity
void compute_precise_delta(
    fe25519* delta,
    const fe25519* z,
    const fe25519* y,
    size_t n
);

// Robust polynomial identity check function
bool robust_polynomial_identity_check(
    const RangeProof* proof,
    const ge25519* V,
    const fe25519* x,
    const fe25519* y,
    const fe25519* z,
    const fe25519* delta,
    const ge25519* g,
    const ge25519* h
);

// Calculate inner product verification point
void calculate_inner_product_point(
    ge25519* P,
    const RangeProof* proof,
    const fe25519* x,
    const fe25519* y,
    const fe25519* z,
    const fe25519* t,
    const PointVector* G,
    const PointVector* H,
    const ge25519* g,
    const ge25519* h,
    size_t n
);

// Verify a range proof
bool range_proof_verify(
    const RangeProof* proof,
    const ge25519* V,       // Value commitment to verify
    size_t n,               // Bit length of range
    const PointVector* G,   // Base points (size n)
    const PointVector* H,   // Base points (size n)
    const ge25519* g,       // Additional base point
    const ge25519* h        // Additional base point
);

void generate_random_scalar(uint8_t* output, size_t len);

// Validate that input is within range
bool validate_range_input(const fe25519* v, size_t n);

// Generate a range proof
void generate_range_proof(
    RangeProof* proof,
    const fe25519* v,
    const fe25519* gamma,
    size_t n,
    const PointVector* G,
    const PointVector* H,
    const ge25519* g,
    const ge25519* h
);

#endif // BULLETPROOF_RANGE_PROOF_H
