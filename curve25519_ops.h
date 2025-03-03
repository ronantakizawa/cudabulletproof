#ifndef CURVE25519_OPS_H
#define CURVE25519_OPS_H

#include <stdint.h>
#include <string.h>

// Field size - 2^255 - 19
#define P25519 0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFED

// Curve constants
#define CURVE25519_A 486662
#define CURVE25519_D 0x52036CEE2B6FFE738CC740797779E89800700A4D4141D8AB75EB4DCA135978A3

// Field element representation for Curve25519
typedef struct {
    uint64_t limbs[4];  // 4 * 64 bit = 256 bits to represent field elements (little-endian)
} fe25519;

// Point representation for Curve25519 in extended coordinates
typedef struct {
    fe25519 X;
    fe25519 Y;
    fe25519 Z;
    fe25519 T;  // X*Y/Z
} ge25519;

// Point in compressed format (for storage/transmission)
typedef struct {
    uint8_t bytes[32];  // Y with sign bit for X in the top bit
} ge25519_compressed;

// Initialize field element from 32-byte array
void fe25519_frombytes(fe25519 *r, const uint8_t *bytes);

// Convert field element to 32-byte array
void fe25519_tobytes(uint8_t *bytes, const fe25519 *h);

// Set field element to 0
void fe25519_0(fe25519 *h);

// Set field element to 1
void fe25519_1(fe25519 *h);

// Copy field element: h = f
void fe25519_copy(fe25519 *h, const fe25519 *f);

// Constant-time conditional swap of field elements
void fe25519_cswap(fe25519 *f, fe25519 *g, uint8_t b);

// Field element addition: h = f + g mod P25519
void fe25519_add(fe25519 *h, const fe25519 *f, const fe25519 *g);

// Field element subtraction: h = f - g mod P25519
void fe25519_sub(fe25519 *h, const fe25519 *f, const fe25519 *g);

// Field element multiplication: h = f * g mod P25519
void fe25519_mul(fe25519 *h, const fe25519 *f, const fe25519 *g);

// Field element squaring: h = f^2 mod P25519
void fe25519_sq(fe25519 *h, const fe25519 *f);

// Field element inversion: h = 1/f mod P25519
void fe25519_invert(fe25519 *h, const fe25519 *f);

// Field element negation: h = -f mod P25519
void fe25519_neg(fe25519 *h, const fe25519 *f);

// Field element power by 2^252 - 3: h = f^(2^252 - 3) mod P25519
// Used in square root computation
void fe25519_pow2523(fe25519 *h, const fe25519 *f);

// Point operations

// Initialize point to identity (neutral element)
void ge25519_0(ge25519 *h);

// Check if point is on curve
int ge25519_is_on_curve(const ge25519 *p);

// Check if point is the identity element
int ge25519_is_identity(const ge25519 *p);

// Point doubling: r = 2*p
void ge25519_double(ge25519 *r, const ge25519 *p);

// Point addition: r = p + q
void ge25519_add(ge25519 *r, const ge25519 *p, const ge25519 *q);

// Point subtraction: r = p - q
void ge25519_sub(ge25519 *r, const ge25519 *p, const ge25519 *q);

// Scalar multiplication: r = scalar * p
void ge25519_scalarmult(ge25519 *r, const uint8_t *scalar, const ge25519 *p);

// Fixed-base scalar multiplication: r = scalar * base
void ge25519_scalarmult_base(ge25519 *r, const uint8_t *scalar);

// Convert point to compressed format
void ge25519_pack(ge25519_compressed *r, const ge25519 *p);

// Convert point from compressed format
int ge25519_unpack(ge25519 *r, const ge25519_compressed *p);

// Copy a point: h = f
void ge25519_copy(ge25519 *h, const ge25519 *f);

// Normalize a point's coordinates to Z=1
void ge25519_normalize(ge25519 *p);

// Device (CUDA) versions of key operations
#ifdef __CUDACC__
__device__ void device_fe25519_add(fe25519 *h, const fe25519 *f, const fe25519 *g);
__device__ void device_fe25519_sub(fe25519 *h, const fe25519 *f, const fe25519 *g);
__device__ void device_fe25519_mul(fe25519 *h, const fe25519 *f, const fe25519 *g);
__device__ void device_fe25519_frombytes(fe25519 *h, const uint8_t *bytes);
__device__ void device_ge25519_add(ge25519 *r, const ge25519 *p, const ge25519 *q);
__device__ void device_ge25519_scalarmult(ge25519 *r, const uint8_t *scalar, const ge25519 *p);
__device__ void device_ge25519_copy(ge25519 *h, const ge25519 *f);
#endif

#endif // CURVE25519_OPS_H
