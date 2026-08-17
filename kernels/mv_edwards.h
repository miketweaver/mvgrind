#ifndef MV_EDWARDS_H
#define MV_EDWARDS_H

typedef struct {
    fe X, Y, Z, T;
} mv_ge;

typedef struct {
    fe Yp, Ym, Z, T2;
} mv_ge_cached;

/* Carry pass. Monocypher's fe_add/fe_sub don't propagate carries, and fe_mul
 * overflows silently (only after ~93 doublings!) if a value chains 2+ adds/subs
 * before a multiply. Never remove a MV_FE_REDUCE without checking that bound. */
#define MV_FE_REDUCE(x) fe_mul_small((x), (x), 1)

static void mv_fe_d2(fe out)
{
    const fe D2 = {-21827239, -5839606, -30745221, 13898782, 229458, 15978800, -12551817, -6495438, 29715968, 9444199};
    fe_copy(out, D2);
}

static void mv_ge_identity(mv_ge *p)
{
    fe_0(p->X);
    fe_1(p->Y);
    fe_1(p->Z);
    fe_0(p->T);
}

static void mv_ge_base(mv_ge *b)
{
    const u8 bx[32] = {0x1a, 0xd5, 0x25, 0x8f, 0x60, 0x2d, 0x56, 0xc9, 0xb2, 0xa7, 0x25, 0x95, 0x60, 0xc7, 0x2c, 0x69,
                       0x5c, 0xdc, 0xd6, 0xfd, 0x31, 0xe2, 0xa4, 0xc0, 0xfe, 0x53, 0x6e, 0xcd, 0xd3, 0x36, 0x69, 0x21};
    const u8 by[32] = {0x58, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
                       0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66};
    fe_frombytes(b->X, bx);
    fe_frombytes(b->Y, by);
    fe_1(b->Z);
    fe_mul(b->T, b->X, b->Y);
}

static void mv_ge_cache(mv_ge_cached *c, const mv_ge *p)
{
    fe d2;
    fe_add(c->Yp, p->Y, p->X);
    fe_sub(c->Ym, p->Y, p->X);
    fe_copy(c->Z, p->Z);
    mv_fe_d2(d2);
    fe_mul(c->T2, p->T, d2);
}

static void mv_ge_add(mv_ge *r, const mv_ge *p, const mv_ge_cached *q)
{
    fe A, B, C, D, t1, t2;

    fe_sub(A, p->Y, p->X);
    fe_mul(A, A, q->Ym);
    fe_add(B, p->Y, p->X);
    fe_mul(B, B, q->Yp);
    fe_mul(C, p->T, q->T2);
    fe_mul(D, p->Z, q->Z);
    fe_add(D, D, D);

    fe_sub(t1, B, A);
    fe_add(B, B, A);
    fe_sub(t2, D, C);
    fe_add(D, D, C);
    MV_FE_REDUCE(t2);
    MV_FE_REDUCE(D);

    fe_mul(r->X, t1, t2);
    fe_mul(r->Y, D, B);
    fe_mul(r->T, t1, B);
    fe_mul(r->Z, t2, D);
}

static void mv_ge_double(mv_ge *r, const mv_ge *p)
{
    fe A, B, C, E, F, G, H;

    fe_sq(A, p->X);
    fe_sq(B, p->Y);
    fe_sq(C, p->Z);
    fe_add(C, C, C);
    fe_add(E, p->X, p->Y);
    fe_sq(E, E);
    fe_sub(E, E, A);
    fe_sub(E, E, B);
    MV_FE_REDUCE(E);
    fe_sub(G, B, A);
    fe_sub(F, G, C);
    MV_FE_REDUCE(F);
    fe_add(H, A, B);
    fe_neg(H, H);

    fe_mul(r->X, E, F);
    fe_mul(r->Y, G, H);
    fe_mul(r->T, E, H);
    fe_mul(r->Z, F, G);
}

static void mv_ge_copy(mv_ge *r, const mv_ge *p)
{
    fe_copy(r->X, p->X);
    fe_copy(r->Y, p->Y);
    fe_copy(r->Z, p->Z);
    fe_copy(r->T, p->T);
}

static void mv_ge_ccopy(mv_ge *r, const mv_ge *p, int b)
{
    fe_ccopy(r->X, p->X, b);
    fe_ccopy(r->Y, p->Y, b);
    fe_ccopy(r->Z, p->Z, b);
    fe_ccopy(r->T, p->T, b);
}

static void mv_ge_scalarmult_base(mv_ge *r, const u8 s[32])
{
    mv_ge b, acc, sum;
    mv_ge_cached bc;

    mv_ge_base(&b);
    mv_ge_cache(&bc, &b);
    mv_ge_identity(&acc);

    for (int i = 255; i >= 0; i--) {
        mv_ge_double(&acc, &acc);
        mv_ge_add(&sum, &acc, &bc);
        mv_ge_ccopy(&acc, &sum, (s[i >> 3] >> (i & 7)) & 1);
    }
    mv_ge_copy(r, &acc);
}

#endif
