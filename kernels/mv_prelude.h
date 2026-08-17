#ifndef MV_PRELUDE_H
#define MV_PRELUDE_H

#if defined(__GNUC__) || defined(__clang__)
#define MV_MAYBE_UNUSED __attribute__((unused))
#else
#define MV_MAYBE_UNUSED
#endif

#ifdef __OPENCL_VERSION__

typedef char i8;
typedef uchar u8;
typedef short i16;
typedef uint u32;
typedef int i32;
typedef long i64;
typedef ulong u64;

#define WIPE_BUFFER(buf) ((void)0)

#else

#include <stddef.h>
#include <stdint.h>

typedef int8_t i8;
typedef uint8_t u8;
typedef int16_t i16;
typedef uint32_t u32;
typedef int32_t i32;
typedef int64_t i64;
typedef uint64_t u64;

static void mv_wipe(void *p, size_t n)
{
    volatile u8 *q = (volatile u8 *)p;
    while (n--)
        *q++ = 0;
}
#define WIPE_BUFFER(buf) mv_wipe((buf), sizeof(buf))

#endif

#define FOR_T(type, i, start, end) for (type i = (start); i < (end); i++)
#define FOR(i, start, end) FOR_T(size_t, i, start, end)
#define COPY(dst, src, size) FOR(_i_, 0, size)(dst)[_i_] = (src)[_i_]
#define ZERO(buf, size) FOR(_i_, 0, size)(buf)[_i_] = 0

static u32 load24_le(const u8 s[3])
{
    return ((u32)s[0] << 0) | ((u32)s[1] << 8) | ((u32)s[2] << 16);
}

static u32 load32_le(const u8 s[4])
{
    return ((u32)s[0] << 0) | ((u32)s[1] << 8) | ((u32)s[2] << 16) | ((u32)s[3] << 24);
}

static void store32_le(u8 out[4], u32 in)
{
    out[0] = (u8)(in);
    out[1] = (u8)(in >> 8);
    out[2] = (u8)(in >> 16);
    out[3] = (u8)(in >> 24);
}

static u64 load64_le(const u8 s[8])
{
    return (u64)load32_le(s) | ((u64)load32_le(s + 4) << 32);
}

static int neq0(u64 diff)
{
    u64 folded = (diff >> 32) | ((u32)diff);
    u64 eq0 = 1 & ((folded - 1) >> 32);
    return (int)eq0 - 1;
}

static u64 x16(const u8 a[16], const u8 b[16])
{
    return (load64_le(a + 0) ^ load64_le(b + 0)) | (load64_le(a + 8) ^ load64_le(b + 8));
}

static u64 x32(const u8 a[32], const u8 b[32])
{
    return x16(a, b) | x16(a + 16, b + 16);
}

static int crypto_verify32(const u8 a[32], const u8 b[32])
{
    return neq0(x32(a, b));
}

static u32 mv_crc32_core(const u8 *buf, u32 len)
{
    u32 crc = 0xFFFFFFFFu;
    for (u32 i = 0; i < len; i++) {
        crc ^= (u32)buf[i];
        for (int k = 0; k < 8; k++) {
            u32 mask = (u32)(-(i32)(crc & 1u));
            crc = (crc >> 1) ^ (0xEDB88320u & mask);
        }
    }
    return ~crc;
}

/* Only ever emit clamped secrets (s = 2^254 + 8k): the firmware derives its
 * XEdDSA signing key from a clamped copy, so an unclamped key produces a node
 * whose signatures don't verify against its own public key. */
static void mv_clamp(u8 s[32])
{
    s[0] &= 0xF8;
    s[31] &= 0x7F;
    s[31] |= 0x40;
}

static int mv_is_clamped(const u8 s[32])
{
    return (s[0] & 0x07) == 0 && (s[31] & 0x80) == 0 && (s[31] & 0x40) == 0x40;
}

MV_MAYBE_UNUSED static void mv_step8(u8 s[32])
{
    u32 carry = 8;
    for (int i = 0; i < 32 && carry; i++) {
        u32 v = (u32)s[i] + carry;
        s[i] = (u8)v;
        carry = v >> 8;
    }
}

static void mv_add8k(u8 s[32], u64 k)
{
    u64 lo = k << 3;
    u64 hi = k >> 61;
    u32 carry = 0;
    for (int i = 0; i < 8; i++) {
        u32 v = (u32)s[i] + (u32)((lo >> (8 * i)) & 0xFFu) + carry;
        s[i] = (u8)v;
        carry = v >> 8;
    }
    for (int i = 8; i < 16; i++) {
        u32 v = (u32)s[i] + (u32)((hi >> (8 * (i - 8))) & 0xFFu) + carry;
        s[i] = (u8)v;
        carry = v >> 8;
    }
    for (int i = 16; i < 32 && carry; i++) {
        u32 v = (u32)s[i] + carry;
        s[i] = (u8)v;
        carry = v >> 8;
    }
}

#endif
