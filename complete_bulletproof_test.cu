
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <openssl/rand.h>
#include <openssl/sha.h>  // Added for SHA256 functions
#include <openssl/crypto.h>  // Added for OPENSSL_init_crypto

#include "curve25519_ops.h"
#include "bulletproof_vectors.h"
#include "bulletproof_range_proof.h"
#include "cuda_bulletproof.h"  // Added to include CUDA function declarations

// External helper logging functions (declared in bulletproof_range_proof.cu)
extern void print_field_element(const char* label, const fe25519* f);
extern void print_point(const char* label, const ge25519* p);

// Forward declarations of our enhanced implementations
void generate_random_scalar(uint8_t* output, size_t len);
void generate_challenge(uint8_t* output, const void* data, size_t data_len, const char* domain_sep);
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

// Generate a deterministic set of base points (in practice, these should be generated in a trusted setup)
void generate_deterministic_base_points(PointVector* points, size_t n, uint8_t seed[32]) {
    for (size_t i = 0; i < n; i++) {
        // Use seed + index to deterministically generate points
        uint8_t hash_input[36];
        memcpy(hash_input, seed, 32);
        hash_input[32] = (i >> 24) & 0xFF;
        hash_input[33] = (i >> 16) & 0xFF;
        hash_input[34] = (i >> 8) & 0xFF;
        hash_input[35] = i & 0xFF;

        // Hash the input to get deterministic bytes
        uint8_t point_bytes[64];
        SHA256_CTX sha_ctx;
        SHA256_Init(&sha_ctx);
        SHA256_Update(&sha_ctx, hash_input, sizeof(hash_input));
        SHA256_Final(point_bytes, &sha_ctx);

        // Use another hash for the Y coordinate
        SHA256_Init(&sha_ctx);
        SHA256_Update(&sha_ctx, point_bytes, 32);
        SHA256_Final(point_bytes + 32, &sha_ctx);

        // Set coordinates
        fe25519_frombytes(&points->elements[i].X, point_bytes);
        fe25519_frombytes(&points->elements[i].Y, point_bytes + 32);

        // Set Z to 1 and compute T = X*Y (for proper curve point)
        fe25519_1(&points->elements[i].Z);
        fe25519_mul(&points->elements[i].T, &points->elements[i].X, &points->elements[i].Y);
    }
}

int main() {
    // Initialize OpenSSL
    OPENSSL_init_crypto(0, NULL);

    // Set bit length for the range proof
    int range_bits = 16;

    printf("Creating a complete Bulletproof range proof with %d bits (CUDA-accelerated)\n", range_bits);

    // Create and initialize base points deterministically
    PointVector G, H;
    point_vector_init(&G, range_bits);
    point_vector_init(&H, range_bits);

    uint8_t G_seed[32] = {0x01};
    uint8_t H_seed[32] = {0x02};
    generate_deterministic_base_points(&G, range_bits, G_seed);
    generate_deterministic_base_points(&H, range_bits, H_seed);

    // Create base points g and h deterministically
    ge25519 g, h;
    uint8_t g_seed[32] = {0x03};
    uint8_t h_seed[32] = {0x04};

    ge25519_0(&g);
    ge25519_0(&h);

    uint8_t g_point_bytes[64], h_point_bytes[64];
    SHA256_CTX sha_ctx;
    SHA256_Init(&sha_ctx);
    SHA256_Update(&sha_ctx, g_seed, 32);
    SHA256_Final(g_point_bytes, &sha_ctx);

    SHA256_Init(&sha_ctx);
    SHA256_Update(&sha_ctx, h_seed, 32);
    SHA256_Final(h_point_bytes, &sha_ctx);

    fe25519_frombytes(&g.X, g_point_bytes);
    fe25519_frombytes(&h.X, h_point_bytes);
    fe25519_1(&g.Y);
    fe25519_1(&h.Y);
    fe25519_1(&g.Z);
    fe25519_1(&h.Z);
    fe25519_mul(&g.T, &g.X, &g.Y);
    fe25519_mul(&h.T, &h.X, &h.Y);

    // Create a value in the range [0, 2^range_bits)
    fe25519 value;
    uint8_t value_bytes[32] = {0};

    // Set a specific value (e.g., 42) within range
    value_bytes[0] = 42;  // Value = 42 (well within range of 2^16)
    fe25519_frombytes(&value, value_bytes);

    // Print the value being tested
    printf("\nTesting value: %d (in range: 0 to %d)\n", value_bytes[0], (1 << range_bits) - 1);

    // Create a random blinding factor
    fe25519 blinding;
    uint8_t blinding_bytes[32];
    generate_random_scalar(blinding_bytes, 32);
    fe25519_frombytes(&blinding, blinding_bytes);

    printf("\nBlinding factor (first 8 bytes): ");
    for (int i = 0; i < 8; i++) {
        printf("%02x", blinding_bytes[i]);
    }
    printf("...\n");

    // Create a value commitment
    ge25519 V;
    pedersen_commit(&V, &value, &blinding, &g, &h);

    printf("Value commitment generated.\n");
    print_point("V", &V);

    // Generate a complete range proof
    RangeProof proof;
    printf("\nGenerating complete range proof...\n");
    generate_range_proof(&proof, &value, &blinding, range_bits, &G, &H, &g, &h);
    printf("Proof generation complete.\n");

    // Verify the range proof using CUDA-accelerated verification
    printf("\n========= STARTING CUDA VERIFICATION =========\n");
    // Benchmark the time difference
    clock_t start_time = clock();

    // Use CUDA-optimized verification
    bool verified = cuda_range_proof_verify(&proof, &V, range_bits, &G, &H, &g, &h);

    clock_t end_time = clock();
    double cuda_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("CUDA Verification Time: %.6f seconds\n", cuda_time);

    printf("========= CUDA VERIFICATION COMPLETE =========\n");

    // For comparison, also run CPU verification
    printf("\n========= STARTING CPU VERIFICATION =========\n");
    start_time = clock();

    // Use CPU verification
    bool cpu_verified = range_proof_verify(&proof, &V, range_bits, &G, &H, &g, &h);

    end_time = clock();
    double cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("CPU Verification Time: %.6f seconds\n", cpu_time);
    printf("========= CPU VERIFICATION COMPLETE =========\n");

    // Print speedup
    if (cpu_time > 0) {
        printf("CUDA Speedup: %.2fx\n", cpu_time / cuda_time);
    }

    // Print result
    printf("\nCUDA Verification result: %s\n", verified ? "SUCCESS" : "FAILED");
    printf("CPU Verification result: %s\n", cpu_verified ? "SUCCESS" : "FAILED");

    if (verified) {
        printf("Successfully verified that the value is in range [0, 2^%d).\n", range_bits);
    } else {
        printf("Verification failed. This could indicate an implementation issue or an invalid value.\n");
        printf("Possible issues to check:\n");
        printf("1. Challenge generation consistency\n");
        printf("2. Point and field element arithmetic\n");
        printf("3. Polynomial coefficient computation\n");
        printf("4. Inner product computation\n");
    }

    // Try with a value outside the range to confirm negative case
    printf("\nTesting with a value outside the range...\n");

    // Create a value outside the range [0, 2^range_bits)
    fe25519 large_value;
    uint8_t large_value_bytes[32] = {0};

    // Set value to 2^range_bits (just outside valid range)
    // For 16-bit range, this would be 65536 (0x10000)
    // We need to handle this explicitly for clarity
    if (range_bits == 16) {
        // little-endian representation of 65536 (0x10000)
        large_value_bytes[0] = 0x00;
        large_value_bytes[1] = 0x00;
        large_value_bytes[2] = 0x01; // This is the 3rd byte (for bit 16)
        large_value_bytes[3] = 0x00;
    } else {
        // Generic calculation for other range_bits values
        large_value_bytes[range_bits/8] |= (1 << (range_bits % 8));
    }

    // Print for debugging
    printf("Testing value: %d (outside range: 0 to %d)\n", 1 << range_bits, (1 << range_bits) - 1);
    printf("Value bytes: ");
    for (int i = 0; i < 4; i++) {
        printf("%02x ", large_value_bytes[i]);
    }
    printf("...\n");

    fe25519_frombytes(&large_value, large_value_bytes);

    // Create a new blinding factor
    fe25519 large_blinding;
    uint8_t large_blinding_bytes[32];
    generate_random_scalar(large_blinding_bytes, 32);
    fe25519_frombytes(&large_blinding, large_blinding_bytes);

    printf("\nBlinding factor for large value (first 8 bytes): ");
    for (int i = 0; i < 8; i++) {
        printf("%02x", large_blinding_bytes[i]);
    }
    printf("...\n");

    // Create commitment
    ge25519 large_V;
    pedersen_commit(&large_V, &large_value, &large_blinding, &g, &h);

    // Generate a range proof (which should fail or be invalid)
    RangeProof large_proof;
    printf("Generating range proof for out-of-range value...\n");
    generate_range_proof(&large_proof, &large_value, &large_blinding, range_bits, &G, &H, &g, &h);

    // Verify using CUDA (should fail)
    printf("Verifying range proof for out-of-range value with CUDA...\n");
    bool large_verified = cuda_range_proof_verify(&large_proof, &large_V, range_bits, &G, &H, &g, &h);

    printf("CUDA Verification result for out-of-range value: %s\n", large_verified ? "SUCCESS (INCORRECT!)" : "FAILED (CORRECT)");

    if (!large_verified) {
        printf("Correctly rejected proof for value outside range, as expected.\n");
    } else {
        printf("Warning: Successfully verified a value outside the range! Implementation issue detected.\n");
    }

    // CUDA-optimized benchmark for field operations
    printf("\n========= CUDA FIELD OPERATIONS BENCHMARK =========\n");
    size_t batch_size = 10000;

    // Create test data
    fe25519* a = (fe25519*)malloc(batch_size * sizeof(fe25519));
    fe25519* b = (fe25519*)malloc(batch_size * sizeof(fe25519));
    fe25519* results = (fe25519*)malloc(batch_size * sizeof(fe25519));

    // Initialize with random data
    for (size_t i = 0; i < batch_size; i++) {
        uint8_t a_bytes[32], b_bytes[32];
        generate_random_scalar(a_bytes, 32);
        generate_random_scalar(b_bytes, 32);
        fe25519_frombytes(&a[i], a_bytes);
        fe25519_frombytes(&b[i], b_bytes);
    }

    // Benchmark CUDA field operations
    printf("Benchmarking batch field operations with %zu elements...\n", batch_size);

    // Addition
    start_time = clock();
    cuda_batch_field_add(results, a, b, batch_size);
    end_time = clock();
    printf("CUDA field addition: %.6f seconds\n", ((double)(end_time - start_time)) / CLOCKS_PER_SEC);

    // Multiplication
    start_time = clock();
    cuda_batch_field_mul(results, a, b, batch_size);
    end_time = clock();
    printf("CUDA field multiplication: %.6f seconds\n", ((double)(end_time - start_time)) / CLOCKS_PER_SEC);

    // Squaring
    start_time = clock();
    cuda_batch_field_square(results, a, batch_size);
    end_time = clock();
    printf("CUDA field squaring: %.6f seconds\n", ((double)(end_time - start_time)) / CLOCKS_PER_SEC);

    // Clean up resources
    free(a);
    free(b);
    free(results);

    // Clean up resources from range proof tests
    point_vector_free(&G);
    point_vector_free(&H);
    range_proof_free(&proof);
    range_proof_free(&large_proof);

    return 0;
}
