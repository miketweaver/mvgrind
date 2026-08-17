__kernel void mv_grind(__constant const uchar *base_scalar, const ulong k_base, const uint iters, __constant const uint *targets,
                       const uint n_targets, const uint mask, volatile __global uint *out_count, __global ulong *out_k,
                       __global uint *out_id, const uint out_cap)
{
    const ulong gid = (ulong)get_global_id(0);
    ulong k = k_base + gid * (ulong)iters;

    u8 s[32];
    for (int i = 0; i < 32; i++)
        s[i] = base_scalar[i];
    mv_add8k(s, k);

    for (uint j = 0; j < iters; j++) {
        u8 pk[32];
        crypto_x25519_public_key(pk, s);
        const u32 id = mv_crc32_core(pk, 32);
        const u32 v = id & mask;

        for (uint t = 0; t < n_targets; t++) {
            if (targets[t] == v) {
                const uint slot = atomic_inc(out_count);
                if (slot < out_cap) {
                    out_k[slot] = k + (ulong)j;
                    out_id[slot] = id;
                }
                break;
            }
        }
        mv_step8(s);
    }
}

static void mv_load_ge(mv_ge *p, __global const int *src)
{
    for (int i = 0; i < 10; i++) {
        p->X[i] = src[i];
        p->Y[i] = src[10 + i];
        p->Z[i] = src[20 + i];
        p->T[i] = src[30 + i];
    }
}

static void mv_load_cached(mv_ge_cached *c, __global const int *src)
{
    for (int i = 0; i < 10; i++) {
        c->Yp[i] = src[i];
        c->Ym[i] = src[10 + i];
        c->Z[i] = src[20 + i];
        c->T2[i] = src[30 + i];
    }
}

/* Do NOT add #pragma unroll to the MV_BATCH loops below: NVIDIA's OpenCL
 * compiler recurses ~15000 frames deep and segfaults during clBuildProgram. */
#ifndef MV_BATCH
#define MV_BATCH 8
#endif

__kernel void mv_grind_batch(__global const int *p_start, __global const int *jumps, const uint n_jumps, __global const int *step,
                             const ulong k_base, const uint iters, __constant const uint *targets, const uint n_targets,
                             const uint mask, volatile __global uint *out_count, __global ulong *out_k, __global uint *out_id,
                             const uint out_cap)
{
    const uint t = (uint)get_global_id(0);

    mv_ge P;
    mv_load_ge(&P, p_start);

    for (uint j = 0; j < n_jumps; j++) {
        if (t & (1u << j)) {
            mv_ge_cached c;
            mv_load_cached(&c, jumps + (size_t)j * 40);
            mv_ge_add(&P, &P, &c);
        }
    }

    mv_ge_cached S;
    mv_load_cached(&S, step);

    const ulong k = k_base + (ulong)t * (ulong)iters;

    for (uint m = 0; m < iters; m += MV_BATCH) {
        fe num[MV_BATCH], den[MV_BATCH], part[MV_BATCH];
        fe acc;

        for (uint i = 0; i < MV_BATCH; i++) {
            fe_add(num[i], P.Z, P.Y);
            fe_sub(den[i], P.Z, P.Y);
            mv_ge_add(&P, &P, &S);
        }

        fe_1(acc);
        for (uint i = 0; i < MV_BATCH; i++) {
            fe_copy(part[i], acc);
            fe_mul(acc, acc, den[i]);
        }

        fe_invert(acc, acc);

        for (int i = MV_BATCH - 1; i >= 0; i--) {
            fe inv_i;
            fe_mul(inv_i, acc, part[i]);
            fe_mul(acc, acc, den[i]);
            fe_mul(num[i], num[i], inv_i);
        }

        for (uint i = 0; i < MV_BATCH; i++) {
            u8 pk[32];
            fe_tobytes(pk, num[i]);
            const u32 id = mv_crc32_core(pk, 32);
            const u32 v = id & mask;
            for (uint q = 0; q < n_targets; q++) {
                if (targets[q] == v) {
                    const uint slot = atomic_inc(out_count);
                    if (slot < out_cap) {
                        out_k[slot] = k + (ulong)m + (ulong)i;
                        out_id[slot] = id;
                    }
                    break;
                }
            }
        }
    }
}

__kernel void mv_derive(__global const uchar *secrets, __global uchar *pubs, __global uint *ids, const uint n)
{
    const uint i = (uint)get_global_id(0);
    if (i >= n)
        return;

    u8 sk[32], pk[32];
    for (int j = 0; j < 32; j++)
        sk[j] = secrets[i * 32 + j];

    crypto_x25519_public_key(pk, sk);

    for (int j = 0; j < 32; j++)
        pubs[i * 32 + j] = pk[j];
    ids[i] = mv_crc32_core(pk, 32);
}
