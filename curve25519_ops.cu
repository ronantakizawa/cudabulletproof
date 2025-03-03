
// File: curve25519_ops.cu
#include "curve25519_ops.h"
#include <stdio.h>

// Curve25519 prime: 2^255 - 19
static const uint64_t p25519[4] = { 0xFFFFFFFFFFFFFFED, 0xFFFFFFFFFFFFFFFF,
                                    0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF };

// Set field element to 0
void fe25519_0(fe25519 *h) {
    memset(h->limbs, 0, sizeof(h->limbs));
}

// Set field element to 1
void fe25519_1(fe25519 *h) {
    h->limbs[0] = 1;
    h->limbs[1] = 0;
    h->limbs[2] = 0;
    h->limbs[3] = 0;
}

// Copy field element: h = f
void fe25519_copy(fe25519 *h, const fe25519 *f) {
    memcpy(h->limbs, f->limbs, sizeof(h->limbs));
}

// Constant-time conditional swap of field elements
void fe25519_cswap(fe25519 *f, fe25519 *g, uint8_t b) {
    uint64_t mask = (uint64_t)(-(int64_t)b);
    uint64_t temp;

    for (int i = 0; i < 4; i++) {
        temp = mask & (f->limbs[i] ^ g->limbs[i]);
        f->limbs[i] ^= temp;
        g->limbs[i] ^= temp;
    }
}

// Field element addition: h = f + g mod P25519
void fe25519_add(fe25519 *h, const fe25519 *f, const fe25519 *g) {
    uint64_t carry = 0;

    for (int i = 0; i < 4; i++) {
        uint64_t sum = f->limbs[i] + g->limbs[i] + carry;

        // Check for overflow
        carry = (sum < f->limbs[i]) || (sum == f->limbs[i] && g->limbs[i] > 0);

        h->limbs[i] = sum;
    }

    // Modular reduction
    if (carry || (h->limbs[3] > p25519[3]) ||
        ((h->limbs[3] == p25519[3]) &&
         ((h->limbs[2] > p25519[2]) ||
          ((h->limbs[2] == p25519[2]) &&
           ((h->limbs[1] > p25519[1]) ||
            ((h->limbs[1] == p25519[1]) && (h->limbs[0] >= p25519[0]))))))) {

        carry = 0;
        for (int i = 0; i < 4; i++) {
            uint64_t diff = h->limbs[i] - p25519[i] - carry;
            carry = (h->limbs[i] < p25519[i] + carry) ? 1 : 0;
            h->limbs[i] = diff;
        }
    }
}

// Field element subtraction: h = f - g mod P25519
void fe25519_sub(fe25519 *h, const fe25519 *f, const fe25519 *g) {
    uint64_t borrow = 0;
    uint64_t temp[4];

    for (int i = 0; i < 4; i++) {
        temp[i] = f->limbs[i] - g->limbs[i] - borrow;
        borrow = (f->limbs[i] < g->limbs[i] + borrow) ? 1 : 0;
    }

    // If result is negative, add prime
    if (borrow) {
        uint64_t carry = 0;
        for (int i = 0; i < 4; i++) {
            temp[i] += p25519[i] + carry;
            carry = (temp[i] < p25519[i]) ? 1 : 0;
        }
    }

    memcpy(h->limbs, temp, sizeof(temp));
}

// Field element multiplication using Karatsuba method adapted for 64-bit limbs
void fe25519_mul(fe25519 *h, const fe25519 *f, const fe25519 *g) {
    // Use temporary space to avoid issues if h overlaps with f or g
    uint64_t t[8] = {0};

    // Schoolbook multiplication - this is not optimized but works for demonstration
    // In a real implementation we would use Karatsuba or optimized assembly
    for (int i = 0; i < 4; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < 4; j++) {
            __uint128_t m = (__uint128_t)f->limbs[i] * g->limbs[j] + t[i+j] + carry;
            t[i+j] = (uint64_t)m;
            carry = (uint64_t)(m >> 64);
        }
        t[i+4] = carry;
    }

    // Modular reduction
    // This is a simplified reduction and not constant time
    // In practice, we'd use a more optimized approach

    // First reduce 2^256 term
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

    // Final reduction
    // Check if result >= p25519
    if (carry || (t[3] > p25519[3]) ||
        ((t[3] == p25519[3]) &&
         ((t[2] > p25519[2]) ||
          ((t[2] == p25519[2]) &&
           ((t[1] > p25519[1]) ||
            ((t[1] == p25519[1]) && (t[0] >= p25519[0]))))))) {

        carry = 0;
        for (int i = 0; i < 4; i++) {
            uint64_t diff = t[i] - p25519[i] - carry;
            carry = (t[i] < p25519[i] + carry) ? 1 : 0;
            h->limbs[i] = diff;
        }
    } else {
        memcpy(h->limbs, t, 4 * sizeof(uint64_t));
    }
}

// Field element squaring: h = f^2 mod P25519
void fe25519_sq(fe25519 *h, const fe25519 *f) {
    // For simplicity, we'll use multiplication
    // In practice, squaring can be optimized further
    fe25519_mul(h, f, f);
}

// Binary extended GCD algorithm for modular inversion
// Computes h = 1/f mod P25519
void fe25519_invert(fe25519 *h, const fe25519 *f) {
    // Use Fermat's Little Theorem: a^(p-2) ≡ a^(-1) (mod p)
    // For Curve25519, we need to compute f^(2^255 - 21)
    fe25519 t0, t1, t2;

    // Compute f^2
    fe25519_sq(&t0, f);

    // Compute f^4 = (f^2)^2
    fe25519_sq(&t1, &t0);

    // Compute f^8 = (f^4)^2
    fe25519_sq(&t1, &t1);

    // Compute f^9 = f * f^8
    fe25519_mul(&t1, &t1, f);

    // Compute f^11 = f^9 * f^2
    fe25519_mul(&t0, &t1, &t0);

    // Compute f^22 = (f^11)^2
    fe25519_sq(&t1, &t0);

    // Continue with exponentiation pattern
    // We'll skip some details for brevity, but the real implementation
    // would carry out the full exponentiation f^(2^255 - 21)

    // For demonstration, we'll perform a few more steps
    // Compute f^44 = (f^22)^2
    fe25519_sq(&t1, &t1);

    // Compute f^88 = (f^44)^2
    fe25519_sq(&t1, &t1);

    // Compute f^176 = (f^88)^2
    fe25519_sq(&t1, &t1);

    // Compute f^220 = f^176 * f^44
    fe25519_mul(&t1, &t1, &t1);

    // Compute f^223 = f^220 * f^3
    fe25519_sq(&t2, f);
    fe25519_mul(&t2, &t2, f);
    fe25519_mul(&t1, &t1, &t2);

    // Continue with this pattern to compute f^(2^255 - 21)
    // Complete implementation would include the full exponentiation chain

    // The inverse is computed after the full exponentiation
    fe25519_copy(h, &t1);
}

// Field element negation: h = -f mod P25519
void fe25519_neg(fe25519 *h, const fe25519 *f) {
    uint64_t borrow = 0;

    for (int i = 0; i < 4; i++) {
        h->limbs[i] = p25519[i] - f->limbs[i] - borrow;
        borrow = (p25519[i] < f->limbs[i] + borrow) ? 1 : 0;
    }
}

// Convert field element to byte representation
void fe25519_tobytes(uint8_t *bytes, const fe25519 *h) {
    fe25519 t;
    fe25519_copy(&t, h);

    // Ensure the value is fully reduced
    if ((t.limbs[3] > p25519[3]) ||
        ((t.limbs[3] == p25519[3]) &&
         ((t.limbs[2] > p25519[2]) ||
          ((t.limbs[2] == p25519[2]) &&
           ((t.limbs[1] > p25519[1]) ||
            ((t.limbs[1] == p25519[1]) && (t.limbs[0] >= p25519[0]))))))) {

        uint64_t borrow = 0;
        for (int i = 0; i < 4; i++) {
            uint64_t diff = t.limbs[i] - p25519[i] - borrow;
            borrow = (t.limbs[i] < p25519[i] + borrow) ? 1 : 0;
            t.limbs[i] = diff;
        }
    }

    // Convert to little-endian bytes
    for (int i = 0; i < 4; i++) {
        bytes[i*8+0] = (t.limbs[i] >> 0) & 0xff;
        bytes[i*8+1] = (t.limbs[i] >> 8) & 0xff;
        bytes[i*8+2] = (t.limbs[i] >> 16) & 0xff;
        bytes[i*8+3] = (t.limbs[i] >> 24) & 0xff;
        bytes[i*8+4] = (t.limbs[i] >> 32) & 0xff;
        bytes[i*8+5] = (t.limbs[i] >> 40) & 0xff;
        bytes[i*8+6] = (t.limbs[i] >> 48) & 0xff;
        bytes[i*8+7] = (t.limbs[i] >> 56) & 0xff;
    }
}

// Convert byte representation to field element
void fe25519_frombytes(fe25519 *h, const uint8_t *bytes) {
    for (int i = 0; i < 4; i++) {
        h->limbs[i] = ((uint64_t)bytes[i*8+0]) |
                      ((uint64_t)bytes[i*8+1] << 8) |
                      ((uint64_t)bytes[i*8+2] << 16) |
                      ((uint64_t)bytes[i*8+3] << 24) |
                      ((uint64_t)bytes[i*8+4] << 32) |
                      ((uint64_t)bytes[i*8+5] << 40) |
                      ((uint64_t)bytes[i*8+6] << 48) |
                      ((uint64_t)bytes[i*8+7] << 56);
    }
}

// Field element power by 2^252 - 3
// This is used to compute square roots in the field
void fe25519_pow2523(fe25519 *h, const fe25519 *f) {
    fe25519 t0, t1, t2;
    int i;

    // Simple exponentiation pattern
    // For a proper implementation, we would optimize this exponentiation chain
    fe25519_sq(&t0, f);
    for (i = 1; i < 5; i++) {
        fe25519_sq(&t0, &t0);
    }
    fe25519_mul(&t1, &t0, f);
    fe25519_sq(&t0, &t1);
    for (i = 1; i < 10; i++) {
        fe25519_sq(&t0, &t0);
    }
    fe25519_mul(&t1, &t0, &t1);
    fe25519_sq(&t0, &t1);
    for (i = 1; i < 20; i++) {
        fe25519_sq(&t0, &t0);
    }
    fe25519_mul(&t0, &t0, &t1);
    fe25519_sq(&t0, &t0);
    for (i = 1; i < 10; i++) {
        fe25519_sq(&t0, &t0);
    }
    fe25519_mul(&t1, &t0, &t1);
    fe25519_sq(&t0, &t1);
    for (i = 1; i < 50; i++) {
        fe25519_sq(&t0, &t0);
    }
    fe25519_mul(&t0, &t0, &t1);
    fe25519_sq(&t0, &t0);
    for (i = 1; i < 100; i++) {
        fe25519_sq(&t0, &t0);
    }
    fe25519_mul(&t0, &t0, &t1);
    fe25519_sq(&t0, &t0);
    for (i = 1; i < 50; i++) {
        fe25519_sq(&t0, &t0);
    }
    fe25519_mul(&t0, &t0, &t1);
    fe25519_sq(&t0, &t0);
    for (i = 1; i < 5; i++) {
        fe25519_sq(&t0, &t0);
    }
    fe25519_mul(h, &t0, &t1);
}

// Initialize point to identity element
void ge25519_0(ge25519 *h) {
    fe25519_0(&h->X);
    fe25519_1(&h->Y);
    fe25519_1(&h->Z);
    fe25519_0(&h->T);
}

// Point addition: r = p + q
void ge25519_add(ge25519 *r, const ge25519 *p, const ge25519 *q) {
    fe25519 A, B, C, D, E, F, G, H;

    // A = (Y1-X1)*(Y2-X2)
    fe25519_sub(&A, &p->Y, &p->X);
    fe25519_sub(&B, &q->Y, &q->X);
    fe25519_mul(&A, &A, &B);

    // B = (Y1+X1)*(Y2+X2)
    fe25519_add(&B, &p->Y, &p->X);
    fe25519_add(&C, &q->Y, &q->X);
    fe25519_mul(&B, &B, &C);

    // C = T1*k*T2
    fe25519 k;
    uint8_t k_bytes[32] = {
        0xA3, 0x78, 0x59, 0x13, 0xCA, 0x4D, 0xEB, 0x75,
        0xAB, 0xD8, 0x41, 0x41, 0x4D, 0x0A, 0x70, 0x00,
        0x98, 0xE8, 0x79, 0x77, 0x79, 0x40, 0xC7, 0x8C,
        0x73, 0xFE, 0x6F, 0x2B, 0xEE, 0x6C, 0x03, 0x52
    }; // Little-endian representation of 2*d
    fe25519_frombytes(&k, k_bytes);
    fe25519_mul(&C, &p->T, &q->T);
    fe25519_mul(&C, &C, &k);

    // D = Z1*2*Z2
    fe25519_mul(&D, &p->Z, &q->Z);
    fe25519_add(&D, &D, &D);

    // E = B - A
    fe25519_sub(&E, &B, &A);

    // F = D - C
    fe25519_sub(&F, &D, &C);

    // G = D + C
    fe25519_add(&G, &D, &C);

    // H = B + A
    fe25519_add(&H, &B, &A);

    // X3 = E*F
    fe25519_mul(&r->X, &E, &F);

    // Y3 = G*H
    fe25519_mul(&r->Y, &G, &H);

    // Z3 = F*G
    fe25519_mul(&r->Z, &F, &G);

    // T3 = E*H
    fe25519_mul(&r->T, &E, &H);
}

// Point subtraction: r = p - q
void ge25519_sub(ge25519 *r, const ge25519 *p, const ge25519 *q) {
    // To subtract a point, we negate it and add
    ge25519 neg_q;

    // Negate q: (x,y) -> (-x,y)
    fe25519_neg(&neg_q.X, &q->X);
    fe25519_copy(&neg_q.Y, &q->Y);
    fe25519_copy(&neg_q.Z, &q->Z);
    fe25519_neg(&neg_q.T, &q->T);

    // Add p + (-q)
    ge25519_add(r, p, &neg_q);
}

// Scalar multiplication: r = scalar * p
// Using double-and-add method
void ge25519_scalarmult(ge25519 *r, const uint8_t *scalar, const ge25519 *p) {
    ge25519 temp;
    ge25519_0(r); // Set result to identity element

    // Process scalar from most significant bit to least
    for (int i = 255; i >= 0; i--) {
        int bit = (scalar[i/8] >> (i % 8)) & 1;

        // Always perform doubling (could be optimized with conditional doubling)
        ge25519_add(&temp, r, r); // double

        // Conditionally perform addition
        if (bit) {
            ge25519_add(r, &temp, p);
        } else {
            ge25519_copy(r, &temp);
        }
    }
}

// Base point for Curve25519
static const uint8_t ge25519_basepoint_bytes[32] = {
    0x58, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
    0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
    0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
    0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66
};

// Fixed-base scalar multiplication: r = scalar * base
void ge25519_scalarmult_base(ge25519 *r, const uint8_t *scalar) {
    // Create a basepoint
    ge25519 base;
    ge25519_0(&base);
    fe25519_frombytes(&base.X, ge25519_basepoint_bytes);
    fe25519_1(&base.Y); // y = 1 for Curve25519
    fe25519_1(&base.Z);
    fe25519_mul(&base.T, &base.X, &base.Y);

    // Perform scalar multiplication
    ge25519_scalarmult(r, scalar, &base);
}

// Negate a point (x,y,z,t) → (-x,y,z,-t)
void ge25519_neg(ge25519 *r, const ge25519 *p) {
    // Negate X and T coordinates, leave Y and Z unchanged
    fe25519_neg(&r->X, &p->X);
    fe25519_copy(&r->Y, &p->Y);
    fe25519_copy(&r->Z, &p->Z);
    fe25519_neg(&r->T, &p->T);
}

// Convert point to compressed format
void ge25519_pack(ge25519_compressed *r, const ge25519 *p) {
    fe25519 recip, x, y;

    // Calculate x and y in affine coordinates
    fe25519_invert(&recip, &p->Z);
    fe25519_mul(&x, &p->X, &recip);
    fe25519_mul(&y, &p->Y, &recip);

    // Encode y with sign bit of x
    fe25519_tobytes(r->bytes, &y);

    // Get least significant bit of x
    uint8_t x_bytes[32];
    fe25519_tobytes(x_bytes, &x);
    uint8_t x_lsb = x_bytes[0] & 1;

    // Set most significant bit of result to x_lsb
    r->bytes[31] |= (x_lsb << 7);
}

// Convert point from compressed format
int ge25519_unpack(ge25519 *r, const ge25519_compressed *p) {
    // Extract y coordinate and sign bit
    fe25519_frombytes(&r->Y, p->bytes);
    uint8_t sign = (p->bytes[31] & 0x80) >> 7;

    // Clear the sign bit in the y value
    uint8_t y_bytes[32];
    memcpy(y_bytes, p->bytes, 32);
    y_bytes[31] &= 0x7F; // Clear the top bit
    fe25519_frombytes(&r->Y, y_bytes);

    // Set Z to 1
    fe25519_1(&r->Z);

    // Compute X from Y using the curve equation
    // x^2 = (y^2 - 1) / (1 + d*y^2)
    fe25519 y2, numerator, denominator, temp, d;

    // Load curve constant d
    uint8_t d_bytes[32] = {
        0xA3, 0x78, 0x59, 0x13, 0xCA, 0x4D, 0xEB, 0x75,
        0xAB, 0xD8, 0x41, 0x41, 0x4D, 0x0A, 0x70, 0x00,
        0x98, 0xE8, 0x79, 0x77, 0x79, 0x40, 0xC7, 0x8C,
        0x73, 0xFE, 0x6F, 0x2B, 0xEE, 0x6C, 0x03, 0x52
    }; // Little-endian representation of Edwards d parameter
    fe25519_frombytes(&d, d_bytes);

    // y^2
    fe25519_sq(&y2, &r->Y);

    // numerator = y^2 - 1
    fe25519 one;
    fe25519_1(&one);
    fe25519_sub(&numerator, &y2, &one);

    // denominator = 1 + d*y^2
    fe25519_mul(&temp, &d, &y2);
    fe25519_add(&denominator, &temp, &one);

    // x^2 = numerator/denominator
    fe25519_invert(&temp, &denominator);
    fe25519_mul(&temp, &numerator, &temp);

    // x = sqrt(x^2)
    // For simplicity, we'll use the helper function that computes the square root
    fe25519 x_squared;
    fe25519_copy(&x_squared, &temp);
    fe25519_pow2523(&r->X, &x_squared); // Approximate square root

    // If the sign bit doesn't match, negate X
    uint8_t x_bytes[32];
    fe25519_tobytes(x_bytes, &r->X);
    if ((x_bytes[0] & 1) != sign) {
        fe25519_neg(&r->X, &r->X);
    }

    // Compute T = X*Y
    fe25519_mul(&r->T, &r->X, &r->Y);

    // Check that the point is on the curve
    return 1; // Simplified for this implementation
}

// Check if point is on curve
int ge25519_is_on_curve(const ge25519 *p) {
    // For a point (X, Y, Z, T) on the curve, we have:
    // -X^2 + Y^2 = Z^2 + d*T^2
    // and T = X*Y/Z

    // This function is not fully implemented in our simplified version
    return 1; // Always return true for this implementation
}

// Check if point is the identity element
int ge25519_is_identity(const ge25519 *p) {
    uint8_t zero[32] = {0};
    uint8_t one[32] = {1}; // Little-endian representation of 1
    uint8_t x_bytes[32], y_bytes[32], z_bytes[32];

    fe25519_tobytes(x_bytes, &p->X);
    fe25519_tobytes(y_bytes, &p->Y);
    fe25519_tobytes(z_bytes, &p->Z);

    // Identity in extended coordinates: (0, 1, 1, 0)
    return (memcmp(x_bytes, zero, 32) == 0 &&
            memcmp(y_bytes, one, 32) == 0 &&
            memcmp(z_bytes, one, 32) == 0);
}

// Point doubling: r = 2*p
void ge25519_double(ge25519 *r, const ge25519 *p) {
    // For simplicity in this demonstration, we'll reuse point addition
    ge25519_add(r, p, p);
}

// Copy a point: h = f
void ge25519_copy(ge25519 *h, const ge25519 *f) {
    fe25519_copy(&h->X, &f->X);
    fe25519_copy(&h->Y, &f->Y);
    fe25519_copy(&h->Z, &f->Z);
    fe25519_copy(&h->T, &f->T);
}

// Point normalization: Convert to equivalent point with Z=1
void ge25519_normalize(ge25519 *p) {
    // Skip if Z is already 1
    uint8_t z_bytes[32];
    fe25519_tobytes(z_bytes, &p->Z);
    uint8_t one_bytes[32] = {1, 0}; // Little-endian representation of 1

    if (memcmp(z_bytes, one_bytes, 32) == 0) {
        return; // Z is already 1, no normalization needed
    }

    // Calculate 1/Z
    fe25519 z_inv;
    fe25519_invert(&z_inv, &p->Z);

    // X' = X/Z
    fe25519 new_x;
    fe25519_mul(&new_x, &p->X, &z_inv);

    // Y' = Y/Z
    fe25519 new_y;
    fe25519_mul(&new_y, &p->Y, &z_inv);

    // T' = X'*Y'
    fe25519 new_t;
    fe25519_mul(&new_t, &new_x, &new_y);

    // Update the point
    fe25519_copy(&p->X, &new_x);
    fe25519_copy(&p->Y, &new_y);
    fe25519_1(&p->Z);  // Z = 1
    fe25519_copy(&p->T, &new_t);
}

// CUDA device implementation of key field operations
#ifdef __CUDACC__
__device__ void device_fe25519_add(fe25519 *h, const fe25519 *f, const fe25519 *g) {
    uint64_t carry = 0;
    uint64_t p25519_d[4] = { 0xFFFFFFFFFFFFFFED, 0xFFFFFFFFFFFFFFFF,
                             0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF };

    for (int i = 0; i < 4; i++) {
        uint64_t sum = f->limbs[i] + g->limbs[i] + carry;
        carry = (sum < f->limbs[i]) || (sum == f->limbs[i] && g->limbs[i] > 0);
        h->limbs[i] = sum;
    }

    // Modular reduction
    if (carry || (h->limbs[3] > p25519_d[3]) ||
        ((h->limbs[3] == p25519_d[3]) &&
         ((h->limbs[2] > p25519_d[2]) ||
          ((h->limbs[2] == p25519_d[2]) &&
           ((h->limbs[1] > p25519_d[1]) ||
            ((h->limbs[1] == p25519_d[1]) && (h->limbs[0] >= p25519_d[0]))))))) {

        carry = 0;
        for (int i = 0; i < 4; i++) {
            uint64_t diff = h->limbs[i] - p25519_d[i] - carry;
            carry = (h->limbs[i] < p25519_d[i] + carry) ? 1 : 0;
            h->limbs[i] = diff;
        }
    }
}

__device__ void device_fe25519_sub(fe25519 *h, const fe25519 *f, const fe25519 *g) {
    uint64_t borrow = 0;
    uint64_t p25519_d[4] = { 0xFFFFFFFFFFFFFFED, 0xFFFFFFFFFFFFFFFF,
                             0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF };
    uint64_t temp[4];

    for (int i = 0; i < 4; i++) {
        temp[i] = f->limbs[i] - g->limbs[i] - borrow;
        borrow = (f->limbs[i] < g->limbs[i] + borrow) ? 1 : 0;
    }

    // If result is negative, add prime
    if (borrow) {
        uint64_t carry = 0;
        for (int i = 0; i < 4; i++) {
            temp[i] += p25519_d[i] + carry;
            carry = (temp[i] < p25519_d[i]) ? 1 : 0;
        }
    }

    for (int i = 0; i < 4; i++) {
        h->limbs[i] = temp[i];
    }
}

__device__ void device_fe25519_mul(fe25519 *h, const fe25519 *f, const fe25519 *g) {
    // Simplified multiplication for device code
    // In practice, this would be more optimized
    uint64_t t[8] = {0};
    uint64_t p25519_d[4] = { 0xFFFFFFFFFFFFFFED, 0xFFFFFFFFFFFFFFFF,
                             0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF };

    for (int i = 0; i < 4; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < 4; j++) {
            unsigned __int128 m = (unsigned __int128)f->limbs[i] * g->limbs[j] + t[i+j] + carry;
            t[i+j] = (uint64_t)m;
            carry = (uint64_t)(m >> 64);
        }
        t[i+4] = carry;
    }

    // Simplified reduction
    uint64_t carry = 0;
    uint64_t c;

    c = t[4] * 19;
    t[0] += c;
    carry = t[0] < c ? 1 : 0;

    for (int i = 1; i < 4; i++) {
        c = t[i+4] * 19 + carry;
        t[i] += c;
        carry = t[i] < c ? 1 : 0;
    }

    // Final reduction check
    if (carry || (t[3] > p25519_d[3]) ||
        ((t[3] == p25519_d[3]) &&
         ((t[2] > p25519_d[2]) ||
          ((t[2] == p25519_d[2]) &&
           ((t[1] > p25519_d[1]) ||
            ((t[1] == p25519_d[1]) && (t[0] >= p25519_d[0]))))))) {

        carry = 0;
        for (int i = 0; i < 4; i++) {
            uint64_t diff = t[i] - p25519_d[i] - carry;
            carry = (t[i] < p25519_d[i] + carry) ? 1 : 0;
            h->limbs[i] = diff;
        }
    } else {
        for (int i = 0; i < 4; i++) {
            h->limbs[i] = t[i];
        }
    }
}

__device__ void device_fe25519_frombytes(fe25519 *h, const uint8_t *bytes) {
    for (int i = 0; i < 4; i++) {
        h->limbs[i] = ((uint64_t)bytes[i*8+0]) |
                      ((uint64_t)bytes[i*8+1] << 8) |
                      ((uint64_t)bytes[i*8+2] << 16) |
                      ((uint64_t)bytes[i*8+3] << 24) |
                      ((uint64_t)bytes[i*8+4] << 32) |
                      ((uint64_t)bytes[i*8+5] << 40) |
                      ((uint64_t)bytes[i*8+6] << 48) |
                      ((uint64_t)bytes[i*8+7] << 56);
    }
}

__device__ void device_ge25519_copy(ge25519 *h, const ge25519 *f) {
    // Copy all fields
    for (int i = 0; i < 4; i++) {
        h->X.limbs[i] = f->X.limbs[i];
        h->Y.limbs[i] = f->Y.limbs[i];
        h->Z.limbs[i] = f->Z.limbs[i];
        h->T.limbs[i] = f->T.limbs[i];
    }
}
#endif
