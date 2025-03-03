
#ifndef BULLETPROOF_CHALLENGE_H
#define BULLETPROOF_CHALLENGE_H

#include "curve25519_ops.h"
#include <stdlib.h>
#include <string.h>
#include <openssl/sha.h>

// Challenge generation function (declared as extern)
extern void generate_challenge(uint8_t* output, const void* data, size_t data_len, const char* domain_sep);

// Generate y challenge from V, A, S
extern void generate_y_challenge(uint8_t* output, const ge25519* V, const ge25519* A, const ge25519* S);

// Generate z challenge from y challenge
extern void generate_z_challenge(uint8_t* output, const uint8_t* y_challenge);

// Generate x challenge from T1, T2
extern void generate_x_challenge(uint8_t* output, const ge25519* T1, const ge25519* T2);

#endif // BULLETPROOF_CHALLENGE_H
