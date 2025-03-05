
// bulletproof_challenge.cu
#include "bulletproof_challenge.h"

// Deterministic challenge generation for Fiat-Shamir
void generate_challenge(uint8_t* output, const void* data, size_t data_len, const char* domain_sep) {
    SHA256_CTX sha_ctx;
    SHA256_Init(&sha_ctx);

    // Add domain separator
    SHA256_Update(&sha_ctx, domain_sep, strlen(domain_sep));

    // Add data
    SHA256_Update(&sha_ctx, data, data_len);

    // Finalize
    SHA256_Final(output, &sha_ctx);

    // Ensure the scalar is in canonical form for curve25519
    output[31] &= 0x7F;  // Clear high bit
}

// Generate y challenge from V, A, S
void generate_challenge_y(uint8_t* output, const ge25519* V, const ge25519* A, const ge25519* S) {
    uint8_t challenge_data[196]; // V(64) + A(64) + S(64) + domain(4)

    // Copy V
    fe25519_tobytes(challenge_data, &V->X);
    fe25519_tobytes(challenge_data + 32, &V->Y);

    // Copy A
    fe25519_tobytes(challenge_data + 64, &A->X);
    fe25519_tobytes(challenge_data + 96, &A->Y);

    // Copy S
    fe25519_tobytes(challenge_data + 128, &S->X);
    fe25519_tobytes(challenge_data + 160, &S->Y);

    // Add domain separator
    memcpy(challenge_data + 192, "y_ch", 4);

    // Generate challenge
    generate_challenge(output, challenge_data, sizeof(challenge_data), "BulletproofYChal");
}

// Generate z challenge from y challenge
void generate_challenge_z(uint8_t* output, const uint8_t* y_challenge) {
    uint8_t challenge_data[36]; // y(32) + domain(4)

    // Copy y challenge
    memcpy(challenge_data, y_challenge, 32);

    // Add domain separator
    memcpy(challenge_data + 32, "z_ch", 4);

    // Generate challenge
    generate_challenge(output, challenge_data, sizeof(challenge_data), "BulletproofZChal");
}

// Generate x challenge from T1, T2
void generate_challenge_x(uint8_t* output, const ge25519* T1, const ge25519* T2) {
    uint8_t challenge_data[132]; // T1(64) + T2(64) + domain(4)

    // Copy T1
    fe25519_tobytes(challenge_data, &T1->X);
    fe25519_tobytes(challenge_data + 32, &T1->Y);

    // Copy T2
    fe25519_tobytes(challenge_data + 64, &T2->X);
    fe25519_tobytes(challenge_data + 96, &T2->Y);

    // Add domain separator
    memcpy(challenge_data + 128, "xchal", 4);

    // Generate challenge
    generate_challenge(output, challenge_data, sizeof(challenge_data), "BulletproofXChal");
}

// Generate inner product challenge
void generate_challenge_inner_product(uint8_t* output, const uint8_t* transcript_data, size_t transcript_len) {
    // Generate challenge with specific domain separation for inner product
    generate_challenge(output, transcript_data, transcript_len, "BulletproofInnerProduct");
}
