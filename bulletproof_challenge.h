
#ifndef BULLETPROOF_CHALLENGE_H
#define BULLETPROOF_CHALLENGE_H

#include "curve25519_ops.h"
#include <stdlib.h>
#include <string.h>
#include <openssl/sha.h>

/**
 * Base challenge generation function using Fiat-Shamir transform
 *
 * @param output Output buffer for the challenge (32 bytes)
 * @param data Input data to hash
 * @param data_len Length of input data
 * @param domain_sep Domain separation string to prevent cross-protocol attacks
 */
void generate_challenge(
    uint8_t* output,
    const void* data,
    size_t data_len,
    const char* domain_sep
);

/**
 * Generate y challenge from V, A, S for Bulletproof range proof
 *
 * @param output Output buffer for the challenge (32 bytes)
 * @param V Value commitment
 * @param A Polynomial commitment A
 * @param S Polynomial commitment S
 */
void generate_challenge_y(
    uint8_t* output,
    const ge25519* V,
    const ge25519* A,
    const ge25519* S
);

/**
 * Generate z challenge from y challenge for Bulletproof range proof
 *
 * @param output Output buffer for the challenge (32 bytes)
 * @param y_challenge Previous y challenge
 */
void generate_challenge_z(
    uint8_t* output,
    const uint8_t* y_challenge
);

/**
 * Generate x challenge from T1, T2 for Bulletproof range proof
 *
 * @param output Output buffer for the challenge (32 bytes)
 * @param T1 Polynomial commitment T1
 * @param T2 Polynomial commitment T2
 */
void generate_challenge_x(
    uint8_t* output,
    const ge25519* T1,
    const ge25519* T2
);

/**
 * Generate inner product challenge for Bulletproof inner product argument
 *
 * @param output Output buffer for the challenge (32 bytes)
 * @param transcript_data Input transcript data
 * @param transcript_len Length of transcript data
 */
void generate_challenge_inner_product(
    uint8_t* output,
    const uint8_t* transcript_data,
    size_t transcript_len
);

#endif // BULLETPROOF_CHALLENGE_H
