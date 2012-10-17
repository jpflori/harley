#undef ulong
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define ulong unsigned long

#include "gf2x.h"

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "padic.h"
#include "padic_poly.h"
#include "qadic_dense.h"

#ifndef FLINT_CPIMPORT
#define FLINT_CPIMPORT "/home/user/FLINT/flint-2/qadic/CPimport.txt"
#endif

#define DEBUG 0

#define NTL_LENGTH
long polyntl[NTL_LENGTH][5] = {
    {24, 4, 3, 1, 0},
    {240, 8, 5, 3, 0},
    {331, 10, 6, 2, 0},
    {1001, 17, 0, 0, 0},
    {1101, 9, 7, 4, 0},
};

#define MGM_LENGTH 37
long polymgm[MGM_LENGTH][9] = {
    {240, 8, 5, 3, 0, 0, 0, 0, 0},
    {384, 8, 7, 6, 4, 3, 2, 1, 0},
    {444, 9, 7, 5, 4, 3, 2, 1, 0},
    {480, 7, 6, 4, 3, 2, 0, 0, 0},
    {500, 8, 6, 5, 2, 1, 0, 0, 0},
    {555, 9, 8, 7, 6, 5, 4, 3, 0},
    {600, 7, 6, 5, 4, 2, 0, 0, 0},
    {665, 10, 9, 7, 6, 5, 3, 1, 0},
    {666, 8, 6, 5, 3, 2, 0, 0, 0},
    {700, 6, 5, 2, 0, 0, 0, 0, 0},
    {777, 10, 9, 8, 5, 4, 0, 0, 0},
    {800, 9, 7, 1, 0, 0, 0, 0, 0},
    {888, 9, 8, 7, 6, 3, 2, 1, 0},
    {900, 1, 0, 0, 0, 0, 0, 0, 0},
    {999, 9, 8, 5, 4, 1, 0, 0, 0},
    {1001, 5, 3, 1, 0, 0, 0, 0, 0},
    {1050, 12, 8, 7, 4, 3, 0, 0, 0},
    {1101, 9, 7, 4, 0, 0, 0, 0, 0},
    {1150, 10, 8, 3, 2, 1, 0, 0, 0},
    {1200, 10, 6, 5, 3, 2, 0, 0, 0},
    {1250, 11, 9, 8, 6, 2, 0, 0, 0},
    {1301, 11, 10, 1, 0, 0, 0, 0, 0},
    {1400, 8, 3, 1, 0, 0, 0, 0, 0},
    {1501, 6, 5, 2, 0, 0, 0, 0, 0},
    {1665, 9, 8, 7, 6, 2, 0, 0, 0},
    {2202, 10, 8, 6, 3, 2, 0, 0, 0},
    {2222, 10, 9, 7, 6, 4, 2, 1, 0},
    {3303, 12, 9, 6, 4, 3, 0, 0, 0},
    {3333, 10, 9, 8, 4, 3, 2, 1, 0},
    {4404, 12, 10, 6, 5, 2, 0, 0, 0},
    {4444, 11, 10, 9, 5, 2, 0, 0, 0},
    {5555, 13, 7, 6, 4, 1, 0, 0, 0},
    {6666, 13, 12, 5, 2, 1, 0, 0, 0},
    {7777, 7, 5, 2, 0, 0, 0, 0, 0},
    {8888, 14, 13, 10, 4, 1, 0, 0, 0},
    {9999, 13, 11, 8, 5, 4, 3, 2, 0},
    {33333, 13, 11, 10, 8, 7, 4, 2, 0},
};

void fmpz_poly_init_ui(fmpz_poly_t rop, long *op)
{
    long *i = op;
    long *j = op + 1;

    fmpz_poly_init2(rop, *i + 1);
    fmpz_set_ui(rop->coeffs + *i, 1);
    _fmpz_poly_set_length(rop, *i + 1);

    for (; *i != 0; i++)
    {
        fmpz_poly_zero_coeffs(rop, *j + 1, *i);
        fmpz_set_ui(rop->coeffs + *j, 1);
        j++;
    }
}

int fmpz_poly_init_conway_mgm(fmpz_poly_t poly,
                              const fmpz_t p, long d)
{
    long i = 0;

    if (fmpz_cmp_ui(p, 2) != 0)
    {
        printf("Exception (fmpz_poly_init_mgm).  Conway polynomials \n");
        printf("are only available for p equal to 2.\n");
        return 1;
    }

    while (polymgm[i][0] != d && i < MGM_LENGTH)
    {
        i++;
    }

    if (i == MGM_LENGTH)
    {
        printf("Exception (fmpz_poly_init_mgm).  Conway polynomial \n");
        printf("not available for this extension degree.\n");
        return 1;
    }

    fmpz_poly_init_ui(poly, polymgm[i]);

    return 0;
}

void qadic_dense_ctx_modulus(padic_poly_t rop, const qadic_dense_ctx_t ctx)
{
    padic_poly_set(rop, ctx->mod);
}

void qadic_dense_ctx_init(qadic_dense_ctx_t ctx, const padic_poly_t mod,
                    const padic_ctx_t pctx, const char *var)
{
    /* Initialize and copy modulus */
    padic_poly_init(ctx->mod);
    padic_poly_set(ctx->mod, mod);

    /* Complete the initialisation of the context */
    padic_ctx_init(&ctx->pctx, pctx->p, pctx->N, pctx->mode);

    /* Precomputed inverse for fast modular reduction */
    qadic_dense_ctx_init_inv(ctx, &ctx->pctx, padic_poly_degree(mod), pctx->N);

    ctx->var = flint_malloc(strlen(var) + 1);
    strcpy(ctx->var, var);

    return;
}

void qadic_dense_ctx_init_reduce(qadic_dense_ctx_t rop, const qadic_dense_ctx_t op, long N)
{
    const long d = qadic_dense_ctx_degree(op);

    int alloc;
    fmpz_t pow;

    /* Initialisation of the p-adic context */
    padic_ctx_init(&rop->pctx, (&op->pctx)->p, N, (&op->pctx)->mode);
    alloc = _padic_ctx_pow_ui(pow, (&rop->pctx)->N, (&rop->pctx));

    /* Initialize and reduce modulus */
    padic_poly_init2(rop->mod, d + 1);
    _fmpz_vec_scalar_mod_fmpz(rop->mod->coeffs, op->mod->coeffs, d + 1, pow);
    _padic_poly_set_length(rop->mod, d + 1);

    padic_poly_init2(rop->invmod, d - 1);
    _fmpz_vec_scalar_mod_fmpz(rop->invmod->coeffs, op->invmod->coeffs, d - 1, pow);
    _padic_poly_set_length(rop->invmod, d - 1);
    _padic_poly_normalise(rop->invmod);

    rop->var = flint_malloc(strlen(op->var) + 1);
    strcpy(rop->var, op->var);

    if (alloc)
        fmpz_clear(pow);

    return;
}

void qadic_dense_ctx_init_reduce_char_2(qadic_dense_ctx_t rop, const qadic_dense_ctx_t op, long N)
{
    const long d = qadic_dense_ctx_degree(op);

    /* Initialisation of the p-adic context */
    padic_ctx_init(&rop->pctx, (&op->pctx)->p, N, (&op->pctx)->mode);

    /* Initialize and reduce modulus */
    padic_poly_init2(rop->mod, d + 1);
    _fmpz_vec_scalar_fdiv_r_2exp(rop->mod->coeffs, op->mod->coeffs, d + 1, N);
    _padic_poly_set_length(rop->mod, d + 1);

    padic_poly_init2(rop->invmod, d - 1);
    _fmpz_vec_scalar_fdiv_r_2exp(rop->invmod->coeffs, op->invmod->coeffs, d - 1, N);
    _padic_poly_set_length(rop->invmod, d - 1);
    _padic_poly_normalise(rop->invmod);

    rop->var = flint_malloc(strlen(op->var) + 1);
    strcpy(rop->var, op->var);

    return;
}

/* Powers of p should be shared */
void _padic_poly_teichmuller_inc_recursive_char_2(fmpz_poly_t delta, const fmpz_poly_t f0,
                                           const fmpz_poly_t f1, const fmpz_poly_t V,
                                           long d, long m, const long *b)
{
    if (m == 1)
    {
        fmpz_poly_neg(delta, V);
        _fmpz_vec_scalar_fdiv_r_2exp(delta->coeffs, delta->coeffs, delta->length, b[0]);
        _fmpz_poly_normalise(delta);
    }
    else
    {
        long *a, h, i, j, k, l, l0, l1, ld, n;
        fmpz_poly_t f0red, delta0, df0, f1red, delta1, df1, W, Delta;

        /* Number of steps */
        n = FLINT_CLOG2(b[0]) + 1;

        /* Precision needed at each step */
        a = flint_malloc(n * sizeof(long));

        i = 0;
        for (a[i = 0] = b[0]; a[i] > 1; i++)
        {
            a[i + 1] = (a[i] + 1) >> 1;
        }

        /* Degree and length */
        l = d + 1;
        l0 = (l + 1) >> 1;
        l1 = (l >> 1);

        fmpz_poly_init2(f0red, l);
        fmpz_poly_init2(delta0, l0);
        fmpz_poly_init2(df0, l);
        fmpz_poly_init2(f1red, l);
        fmpz_poly_init2(delta1, l1);
        fmpz_poly_init2(df1, l);
        fmpz_poly_init2(W, l);
        fmpz_poly_init2(Delta, l);

        /* Lifting */
        i = n - 1;
        j = m - 1;
        {
            ld = V->length;

            _fmpz_vec_scalar_fdiv_r_2exp(delta->coeffs, V->coeffs, ld, a[i]);
            _fmpz_poly_set_length(delta, ld);
            _fmpz_poly_normalise(delta);
        }
        for (i--; i >= 0; i--)
        {
            h = j;
            j--;
            if (b[j] != a[i])
            {
                j--;
            }

            ld = delta->length;

            /* delta_0, delta_1, and products */
            _fmpz_poly_set_length(delta0, (ld + 1) >> 1);
            _fmpz_poly_set_length(delta1, (ld >> 1));

            for (k = 0; k < ld; k++)
                fmpz_set(((k % 2)?delta1:delta0)->coeffs + (k >> 1), delta->coeffs + k);

            _fmpz_poly_normalise(delta0);
            _fmpz_poly_normalise(delta1);

            _fmpz_vec_scalar_fdiv_r_2exp(f0red->coeffs, f0->coeffs, f0->length, a[i]);
            _fmpz_vec_scalar_fdiv_r_2exp(f1red->coeffs, f1->coeffs, f1->length, a[i]);
            _fmpz_poly_set_length(f0red, f0->length);
            _fmpz_poly_set_length(f1red, f1->length);
            _fmpz_poly_normalise(f0red);
            _fmpz_poly_normalise(f1red);

            if (delta0->length != 0)
            {
                if (delta0->length >= f0red->length)
                    _fmpz_poly_mul(df0->coeffs, delta0->coeffs, delta0->length, f0red->coeffs, f0red->length);
                else
                    _fmpz_poly_mul(df0->coeffs, f0red->coeffs, f0red->length, delta0->coeffs, delta0->length);
                _fmpz_poly_set_length(df0, delta0->length + f0red->length - 1);
                _fmpz_vec_scalar_fdiv_r_2exp(df0->coeffs, df0->coeffs, df0->length, a[i]);
            }
            else
                df0->length = 0;
            if (delta1->length != 0)
            {
                if (delta1->length >= f1red->length)
                    _fmpz_poly_mul(df1->coeffs, delta1->coeffs, delta1->length, f1red->coeffs, f1red->length);
                else
                    _fmpz_poly_mul(df1->coeffs, f1red->coeffs, f1red->length, delta1->coeffs, delta1->length);
                _fmpz_poly_set_length(df1, delta1->length + f1red->length);
                _fmpz_vec_scalar_fdiv_r_2exp(df1->coeffs, df1->coeffs, df1->length, a[i]);
                _fmpz_poly_shift_left(df1->coeffs, df1->coeffs, df1->length - 1, 1);
            }
            else
                df1->length = 0;

            if (df0->length != 0 || df1->length != 0)
            {
                _fmpz_poly_sub(df1->coeffs, df1->coeffs, df1->length, df0->coeffs, df0->length);
                _fmpz_poly_set_length(df1, FLINT_MAX(df0->length, df1->length));
                _fmpz_vec_scalar_fdiv_r_2exp(df1->coeffs, df1->coeffs, df1->length, a[i]-1);
                _fmpz_poly_normalise(df1);
                _fmpz_vec_scalar_mul_2exp(df1->coeffs, df1->coeffs, df1->length, 1);
            }

            /* W */
            _fmpz_poly_add(W->coeffs, delta->coeffs, ld, V->coeffs, V->length);
            _fmpz_poly_set_length(W, FLINT_MAX(ld, V->length));
            _fmpz_vec_scalar_fdiv_r_2exp(W->coeffs, W->coeffs, W->length, a[i]);
            _fmpz_poly_normalise(W);

            _fmpz_poly_add(W->coeffs, W->coeffs, W->length, df1->coeffs, df1->length);
            _fmpz_poly_set_length(W, FLINT_MAX(W->length, df1->length));
            _fmpz_vec_scalar_fdiv_r_2exp(W->coeffs, W->coeffs, W->length, a[i]);
            _fmpz_poly_normalise(W);

            _fmpz_vec_scalar_fdiv_q_2exp(W->coeffs, W->coeffs, W->length, a[i + 1]);

            /* Delta */
            if (b[h] == a[i] - a[i + 1])
                _padic_poly_teichmuller_inc_recursive_char_2(Delta, f0, f1, W, d, m - h, b + h);
            else
                _padic_poly_teichmuller_inc_recursive_char_2(Delta, f0, f1, W, d, m - h - 1, b + h + 1);
            _fmpz_vec_scalar_mul_2exp(Delta->coeffs, Delta->coeffs, Delta->length, a[i + 1]);
            /*_fmpz_vec_scalar_fdiv_r_2exp(Delta->coeffs, Delta->coeffs, Delta->length, a[i]);
              _fmpz_poly_normalise(Delta);*/

            /* update delta */
            _fmpz_poly_add(delta->coeffs, delta->coeffs, ld, Delta->coeffs, Delta->length);
            _fmpz_poly_set_length(delta, FLINT_MAX(ld, Delta->length));
            _fmpz_vec_scalar_fdiv_r_2exp(delta->coeffs, delta->coeffs, delta->length, a[i]);
            _fmpz_poly_normalise(delta);
        }

        fmpz_poly_clear(f0red);
        fmpz_poly_clear(delta0);
        fmpz_poly_clear(df0);
        fmpz_poly_clear(delta1);
        fmpz_poly_clear(df1);
        fmpz_poly_clear(f1red);
        fmpz_poly_clear(W);
        fmpz_poly_clear(Delta);
    }
}

/* Partial lifts could be stored rather than recomputed by reduction later */
void _padic_poly_teichmuller(fmpz_poly_t g, const fmpz_poly_t f,
                             const fmpz_t p, long N)
{
    if (N == 1)
    {
        fmpz_poly_set(g, f);
    }
    /*
      Specialized methods should be used for mod operations
      in char 2
    */
    else if (*(p) == 2L)
    {
        long *a, *b, d, h, i, j, k, l, l0, l1, n, m;
        fmpz_poly_t f0, f0_sqr, f1, f1_sqr, V, delta;

        /* Degree and length */
        d = fmpz_poly_degree(f);
        l = d + 1;
        l0 = (l + 1) >> 1;
        l1 = (l >> 1);

        fmpz_poly_init2(f0, l0);
        fmpz_poly_init2(f0_sqr, l);
        fmpz_poly_init2(f1, l1);
        fmpz_poly_init2(f1_sqr, l);
        fmpz_poly_init2(V, l);
        fmpz_poly_init2(delta, l);

        /* Number of steps */
        n = FLINT_CLOG2(N) + 1;

        /* Precision needed at each step */
        a = flint_malloc(n * sizeof(long));

        i = 0;
        for (a[i = 0] = N; a[i] > 1; i++)
        {
            a[i + 1] = (a[i] + 1) >> 1;
        }

        /* Compute precisions needed later */
        b = flint_malloc((2 * n) * sizeof(long));

        i = 0;
        {
            b[i] = a[i];
        }
        for (; !(b[i] & 1L); i++)
        {
            b[i + 1] = a[i + 1];
        }
        for (; b[i] > 2; i += 2)
        {
            b[i + 1] = (b[i] >> 1) + 1;
            b[i + 2] = b[i] >> 1;
        }
        if (b[i] == 2L)
        {
            b[i + 1] = 1L;
            i++;
        }
        m = i + 1;

        /* Lifting */
        j = m - 1;
        i = n - 1;
        {
            fmpz_poly_set(g, f);
        }
        for (i--; i >= 0; i--)
        {
            h = j;
            j--;
            if (b[j] != a[i])
            {
                j--;
            }

            /* f_0, f_1, and squares */
            _fmpz_poly_set_length(f0, l0);
            _fmpz_poly_set_length(f1, l1);

            for (k = 0; k < l; k++)
                fmpz_set(((k % 2)?f1:f0)->coeffs + (k >> 1), g->coeffs + k);

            _fmpz_poly_normalise(f0);
            _fmpz_poly_normalise(f1);

            _fmpz_poly_sqr(f0_sqr->coeffs, f0->coeffs, f0->length);
            _fmpz_poly_sqr(f1_sqr->coeffs, f1->coeffs, f1->length);

            _fmpz_poly_set_length(f0_sqr, (f0->length << 1) - 1);
            _fmpz_poly_set_length(f1_sqr, (f1->length << 1));

            _fmpz_vec_scalar_fdiv_r_2exp(f0_sqr->coeffs, f0_sqr->coeffs, f0_sqr->length, a[i]);
            _fmpz_vec_scalar_fdiv_r_2exp(f1_sqr->coeffs, f1_sqr->coeffs, f1_sqr->length - 1, a[i]);

            _fmpz_poly_shift_left(f1_sqr->coeffs, f1_sqr->coeffs, f1_sqr->length - 1, 1);

            /* V */
            _fmpz_poly_set_length(V, l);

            _fmpz_vec_set(V->coeffs, g->coeffs, l);

            _fmpz_poly_sub(V->coeffs, V->coeffs, V->length, f0_sqr->coeffs, f0_sqr->length);
            /*_fmpz_vec_scalar_fdiv_r_2exp(V->coeffs, V->coeffs, V->length, a[i]);
              _fmpz_poly_normalise(V);*/

            _fmpz_poly_add(V->coeffs, V->coeffs, V->length, f1_sqr->coeffs, f1_sqr->length);
            _fmpz_poly_set_length(V, FLINT_MAX(V->length, f1_sqr->length));
            /*_fmpz_vec_scalar_fdiv_r_2exp(V->coeffs, V->coeffs, V->length, a[i]);
              _fmpz_poly_normalise(V);*/

            _fmpz_vec_scalar_fdiv_q_2exp(V->coeffs, V->coeffs, V->length, a[i + 1]);
            _fmpz_poly_normalise(V);

            /* delta */
            if (b[h] == a[i] - a[i + 1])
                _padic_poly_teichmuller_inc_recursive_char_2(delta, f0, f1, V, d, m - h, b + h);
            else
                _padic_poly_teichmuller_inc_recursive_char_2(delta, f0, f1, V, d, m - h - 1, b + h + 1);

            /* Use vec method */
            _fmpz_vec_scalar_mul_2exp(delta->coeffs, delta->coeffs, delta->length, a[i + 1]);
            /*_fmpz_vec_scalar_fdiv_r_2exp(delta->coeffs, delta->coeffs, delta->length, a[i]);*/
            _fmpz_poly_normalise(delta);

            /* update g */
            _fmpz_poly_add(g->coeffs, g->coeffs, l, delta->coeffs, delta->length);
            _fmpz_vec_scalar_fdiv_r_2exp(g->coeffs, g->coeffs, l, a[i]);
        }

        if (!fmpz_is_one(fmpz_poly_lead(g)))
        {
            _fmpz_vec_neg(g->coeffs, g->coeffs, l);
            _fmpz_vec_scalar_fdiv_r_2exp(g->coeffs, g->coeffs, l, N);
        }

        fmpz_poly_clear(f0);
        fmpz_poly_clear(f0_sqr);
        fmpz_poly_clear(f1);
        fmpz_poly_clear(f1_sqr);
        fmpz_poly_clear(V);
        fmpz_poly_clear(delta);

        flint_free(a);
        flint_free(b);
    }
}

#define L2W(l) (((l) / (sizeof(unsigned long) * CHAR_BIT)) + (((l) % (sizeof(unsigned long) * CHAR_BIT)) ? 1 : 0))

static void
gf2x_clear(unsigned long *A)
{
    flint_free(A);
}

static void
gf2x_copy(unsigned long *R,
          const unsigned long *A, unsigned int wA, long lA)
{
    unsigned int i;

    for (i = 0; i < wA; i++)
        R[i] = A[i];
}

static void
gf2x_normalise(const unsigned long *A, unsigned int *wA, long *lA)
{
    unsigned long bit;

    while (((*wA) > 0) && !A[(*wA) - 1])
    {
        (*wA)--;
        *lA = (*wA) * (sizeof(unsigned long) * CHAR_BIT);
    }

    if ((*wA) == 0)
        *lA = 0;
    else
    {
        if ((*lA) % (sizeof(unsigned long) * CHAR_BIT) == 0)
            bit = 1UL << ((sizeof(unsigned long) * CHAR_BIT) - 1);
        else
            bit = 1UL << (((*lA) % (sizeof(unsigned long) * CHAR_BIT)) - 1);

        while (!(A[(*wA) - 1] & bit))
        {
            (*lA)--;
            bit >>= 1;
        }
    }
}

/* No aliasing */
static void
gf2x_reverse(unsigned long *A, const unsigned long *B, unsigned int wB, long lB)
{
    long i;
    unsigned int w, wA;
    unsigned long bit, bitrev;

    i = lB;
    w = wB - 1;
    bitrev = 1UL;
    wA = 0;
    A[wA] = 0;

    if (i % (sizeof(unsigned long) * CHAR_BIT))
    {
        bit = 1UL << ((i % (sizeof(unsigned long) * CHAR_BIT)) - 1);

        while (bit != 0)
        {
            if (B[w] & bit)
                A[wA] |= bitrev;

            bit >>= 1;
            bitrev <<= 1;
            i--;
        }
        w--;
    }

    while (i != 0)
    {
        bit = 1UL << ((sizeof(unsigned long) * CHAR_BIT) - 1);

        while (bit != 0)
        {
            if (bitrev == 0)
            {
                bitrev = 1UL;
                wA++;
                A[wA] = 0;
            }

            if (B[w] & bit)
                A[wA] |= bitrev;

            bit >>= 1;
            bitrev <<= 1;
            i--;
        }
        w--;
    }
}

/* Aliasing possible */
static void
gf2x_add(unsigned long *R, const unsigned long *A, unsigned int wA,
         const unsigned long *B, unsigned int wB)
{
    if (wA < wB)
    {
        gf2x_add(R, B, wB, A, wA);
    }
    else
    {
        long i;

        for (i = 0; i < wB; i++)
            R[i] = A[i] ^ B[i];

        if (R != A)
        {
            for (; i < wA; i++)
                R[i] = A[i];
        }
    }
}

/* Aliasing possible */
static void
gf2x_sqr(unsigned long *A, const unsigned long *B, unsigned int wB, long lB)
{
    long i;
    unsigned int w, wA;
    unsigned long bit, bitsqr, C;

    i = lB;
    w = wB - 1;
    wA = L2W(2*lB - 1) - 1;
    if (i % (sizeof(unsigned long) * CHAR_BIT))
    {
        bit = 1UL << ((i % (sizeof(unsigned long) * CHAR_BIT)) - 1);
        bitsqr = 1UL << (((2*i - 1) % (sizeof(unsigned long) * CHAR_BIT)) - 1);

        C = B[w];
        A[wA] = 0;

        while (bit != 0)
        {
            if (bitsqr == 0)
            {
                bitsqr = 1UL << ((sizeof(unsigned long) * CHAR_BIT) - 2);
                wA--;
                A[wA] = 0;
            }

            if (C & bit)
                A[wA] |= bitsqr;

            bit >>= 1;
            i--;

            bitsqr >>= 2;
        }
        w--;
        wA--;
    }
    while (i != 0)
    {
        bit = 1UL << ((sizeof(unsigned long) * CHAR_BIT) - 1);
        bitsqr = 1UL << ((sizeof(unsigned long) * CHAR_BIT) - 2);

        C = B[w];
        A[wA] = 0;

        while (bit != 0)
        {
            if (bitsqr == 0)
            {
                bitsqr = 1UL << ((sizeof(unsigned long) * CHAR_BIT) - 2);
                wA--;
                A[wA] = 0;
            }

            if (C & bit)
                A[wA] |= bitsqr;

            bit >>= 1;
            i--;

            bitsqr >>= 2;
        }
        w--;
        wA--;
    }
}

/* No aliasing */
static __inline__ void
gf2x_inv_series_newton(unsigned long *R, const unsigned long *A, unsigned int wA, long lA,
                       long d)
{
    long i = FLINT_CLOG2(d);
    long t = FLINT_CLOG2(sizeof(unsigned long) * CHAR_BIT);
    long j;
    unsigned int wR;
    long lR;

    R[0] = 1;
    for (j = 1; j < L2W(d); j++)
    {
        R[j] = 0;
    }
    wR = 1;
    lR = 1;

    while (t-- && i--)
    {
        gf2x_sqr(R, R, wR, lR);
        gf2x_mul(R, R, wR, A, 1);
        lR <<= 1;
    }

    if (i > 0)
    {
        long *a;
        unsigned long *W;

        a = flint_malloc((i + 1) * sizeof(long));
        for (a[j = 0] = d; a[j] > (sizeof(unsigned long) * CHAR_BIT); j++)
            a[j + 1] = (a[j] + 1) / 2;

        W = flint_malloc((wA + L2W(a[1])) * sizeof(unsigned long));

        while (i--)
        {
            gf2x_mul(W, R, wR, A, wA);
            gf2x_mul(R + wR, R, wR, W + wR, FLINT_MIN(wR, wA));
            lR = a[i];
            wR = L2W(lR);
        }

        flint_free(a);

        flint_free(W);
    }

    {
        unsigned long mask = (1UL << (d % (sizeof(unsigned long) * CHAR_BIT))) - 1;
        R[wR - 1] &= mask;
    }
}


/* No aliasing */
/* wQ >= L2W(lQ) and wR >= wB + wQ >= wA */
static __inline__ void
gf2x_divrem_newton(unsigned long *Q, unsigned long *R,
                   const unsigned long *A, unsigned int wA, long lA,
                   const unsigned long *B, unsigned int wB, long lB)
{
    if (lA >= lB)
    {
        const long lQ = lA - lB + 1;
        const unsigned int wQ = L2W(lQ);
        unsigned long *Arev, *Brev, *Binv, *Qrev;

        Qrev = flint_malloc((2 * wQ) * sizeof(unsigned long));

        Brev = flint_malloc(wB * sizeof(unsigned long));
        gf2x_reverse(Brev, B, wB, lB);

        Binv = flint_malloc((wB + L2W(lQ) + 1) * sizeof(unsigned long));
        gf2x_inv_series_newton(Binv, Brev, wB, lB, lQ);

        Arev = flint_malloc(wA * sizeof(unsigned long));
        gf2x_reverse(Arev, A, wA, lA);

        gf2x_mul(Qrev, Arev, wQ, Binv, wQ);

        gf2x_reverse(Q, Qrev, wQ, lQ);

        gf2x_mul(R, B, wB, Q, wQ);
        gf2x_add(R, A, wA, R, wA);

        flint_free(Arev);
        flint_free(Brev);
        flint_free(Binv);
        flint_free(Qrev);
    }
    else
    {
        gf2x_copy(R, A, wA, lA);
    }
}

/* No aliasing */
static void
gf2x_xgcd_euclidean(unsigned long *G, unsigned long *U, unsigned long *V,
                    const unsigned long *A, unsigned int wA, long lA,
                    const unsigned long *B, unsigned int wB, long lB)
{
    if (lA < lB)
        gf2x_xgcd_euclidean(G, U, V, B, wB, lB, A, wA, lA);
    else
    {
        long i;
        unsigned long *Q, *R, *C, *D, *T;
        unsigned long *M11, *M12, *M21, *M22, *M31, *M32;
        unsigned int wC, wD, wQ, wR, wM11, wM12, wM21, wM22, wM31, wM32;
        long lC, lD, lQ, lR, lM11, lM12, lM21, lM22, lM31, lM32;

        Q = flint_malloc(wA * sizeof(unsigned long));
        R = flint_malloc((wA + 1) * sizeof(unsigned long));
        C = flint_malloc((wA + 1) * sizeof(unsigned long));
        D = flint_malloc((wA + 1) * sizeof(unsigned long));
        M11 = flint_malloc((wA + 1) * sizeof(unsigned long));
        M12 = flint_malloc((wA + 1) * sizeof(unsigned long));
        M21 = flint_malloc((wA + 1) * sizeof(unsigned long));
        M22 = flint_malloc((wA + 1) * sizeof(unsigned long));
        M31 = flint_malloc((wA + 1) * sizeof(unsigned long));
        M32 = flint_malloc((wA + 1) * sizeof(unsigned long));

        gf2x_copy(C, A, wA, lA);
        wC = wA;
        lC = lA;
        gf2x_copy(D, B, wB, lB);
        wD = wB;
        lD = lB;

        M11[0] = 1;
        wM11 = 1;
        lM11 = 1;
        M22[0] = 1;
        wM22 = 1;
        lM22 = 1;

        wM12 = 0;
        lM12 = 0;
        wM21 = 0;
        lM21 = 0;

        while (lD)
        {
            gf2x_divrem_newton(Q, R, C, wC, lC, D, wD, lD);
            lQ = lC - lD + 1;
            wQ = L2W(lQ);
            lR = lD - 1;
            wR = L2W(lR);
            gf2x_normalise(R, &wR, &lR);

            gf2x_mul(M31, Q, wQ, M21, wM21);
            lM31 = lQ + lM21 - 1;
            wM31 = L2W(lM31);
            gf2x_mul(M32, Q, wQ, M22, wM22);
            lM32 = lQ + lM22 - 1;
            wM32 = L2W(lM32);

            gf2x_add(M31, M31, wM31, M11, wM11);
            T = M11;
            M11 = M21;
            M21 = M31;
            M31 = T;
            lM31 = FLINT_MAX(lM31, lM11);
            wM31 = FLINT_MAX(wM31, wM11);
            lM11 = lM21;
            wM11 = wM21;
            lM21 = lM31;
            wM21 = wM31;

            gf2x_add(M32, M32, wM32, M12, wM12);
            T = M12;
            M12 = M22;
            M22 = M32;
            M32 = T;
            lM32 = FLINT_MAX(lM32, lM12);
            wM32 = FLINT_MAX(wM32, wM12);
            lM12 = lM22;
            wM12 = wM22;
            lM22 = lM32;
            wM22 = wM32;
            /*printf("M11[0] %lX, M12[0] %lX, M21[0] %lX, M22[0] %lX\n", M11[0], M12[0], M21[0], M22[0]);
              printf("w11 %u, w12 %u, w21 %u, w22 %u\n", wM11, wM12, wM21, wM22);*/
            T = C;
            C = D;
            D = R;
            R = T;
            lC = lD;
            wC = wD;
            lD = lR;
            wD = wR;
        }

        for (i = 0; i < wB; i++)
            G[i] = 0;

        for (i = 0; i < wA; i++)
        {
            U[i] = 0;
            V[i] = 0;
        }

        gf2x_copy(G, C, wC, lC);
        gf2x_copy(U, M11, wM11, lM11);
        gf2x_copy(V, M12, wM12, lM12);

        flint_free(C);
        flint_free(D);
        flint_free(Q);
        flint_free(R);
        flint_free(M11);
        flint_free(M12);
        flint_free(M21);
        flint_free(M22);
        flint_free(M31);
        flint_free(M32);

        /*printf("G[0] %lX U[0] %lX V[0] %lX\n", G[0], U[0], V[0]);*/
    }
}

/* No aliasing */
static void
gf2x_xgcd(unsigned long *G, unsigned long *U, unsigned long *V,
          const unsigned long *A, unsigned int wA, long lA,
          const unsigned long *B, unsigned int wB, long lB)
{
    if (lA < lB)
        gf2x_xgcd(G, U, V, B, wB, lB, A, wA, lA);
    else
        gf2x_xgcd_euclidean(G, U, V, A, wA, lA, B, wB, lB);

}


/* No aliasing */
/* Suppose lA < lB */
static void
gf2x_invmod(unsigned long *R, const unsigned long *A, unsigned int wA, long lA,
            const unsigned long *B, unsigned int wB, long lB)
{
    unsigned long *S, *T;
    S = flint_malloc(wB * sizeof(unsigned long));
    T = flint_malloc(wB * sizeof(unsigned long));

    gf2x_xgcd(S, T, R, B, wB, lB, A, wA, lA);

    flint_free(S);
    flint_free(T);
}

/* No aliasing */
static __inline__ void
gf2x_reduce(unsigned long *R, const unsigned long *A, unsigned int wA, long lA,
            const unsigned long *B, unsigned int wB, long lB,
            const unsigned long *Binv, unsigned int wBinv, long lBinv)
{
    if (lA >= lB)
    {
        const long lQ = lA - lB + 1;
        const unsigned int wQ = L2W(lQ);
        unsigned long *Arev, *Q, *Qrev;

        Q = flint_malloc(wQ * sizeof(unsigned long));
        Qrev = flint_malloc((2 * wQ) * sizeof(unsigned long));

        Arev = flint_malloc(wA * sizeof(unsigned long));
        gf2x_reverse(Arev, A, wA, lA);

        gf2x_mul(Qrev, Arev, FLINT_MIN(wA, wQ), Binv, FLINT_MIN(wBinv, wQ));
        gf2x_reverse(Q, Qrev, wQ, lQ);

        gf2x_mul(R, B, wB, Q, wQ);
        gf2x_add(R, A, wA, R, wB);

        flint_free(Arev);
        flint_free(Q);
        flint_free(Qrev);
    }
    else
    {
        gf2x_copy(R, A, wA, lA);
    }
}

/* Aliasing possible */
static void
_qadic_dense_sqrt_gf2x(unsigned long *A, unsigned int *wA, long *lA,
                        const unsigned long *B, unsigned int wB, long lB,
                        const unsigned long *sqrt, unsigned int wsqrt, long lsqrt,
                        const unsigned long *mod, unsigned int wmod, long lmod,
                        const unsigned long *invmod, unsigned int winvmod, long linvmod)
{
    long k;
    unsigned long *s;
    unsigned int ws, w;
    long ls;
    unsigned long bit, bitsqrt, t;

    s =  (unsigned long *) flint_malloc((2 * wmod) * sizeof(unsigned long));

    {
        bit = 1UL;
        bitsqrt = 1UL;

        w = 0;
        t = B[w];

        *wA = 0;
        A[*wA] = 0;

        ws = 0;
        s[ws] = 0;

        for (k = 0; k < lB; k++)
        {
            if (bit == 0)
            {
                bit = 1UL;
                w++;
                t = B[w];
            }

            if (bitsqrt == 0)
            {
                bitsqrt = 1UL;
                (*wA)++;
                A[*wA] = 0;
                ws++;
                s[ws] = 0;
            }

            if (k & 1)
            {
                if (t & bit)
                    s[ws] |= bitsqrt;
                bitsqrt <<= 1;
            }
            else
            {
                if (t & bit)
                    A[*wA] |= bitsqrt;
            }

            bit <<=1;
        }

        *lA = (lB + 1) >> 1;
        (*wA)++;
        gf2x_normalise(A, wA, lA);

        ls = lB >> 1;
        ws++;
        gf2x_normalise(s, &ws, &ls);
    }
    {
        gf2x_mul(s, s, ws, sqrt, wsqrt);
        ls = ls + lsqrt - 1;
        ws = L2W(ls);

        gf2x_add(s, s, ws, A, *wA);

        gf2x_reduce(A, s, ws, ls, mod, wmod, lmod, invmod, winvmod, linvmod);
        *lA = FLINT_MIN(lmod - 1, ls);
        *wA = L2W(*lA);
        gf2x_normalise(A, wA, lA);
    }

    flint_free(s);
}

/* Aliasing possible */
static
void _qadic_dense_sqr_gf2x(unsigned long *A, unsigned int *wA, long *lA,
                           const unsigned long *B, unsigned int wB, long lB,
                           const unsigned long *mod, unsigned int wmod, long lmod,
                           const unsigned long *invmod, unsigned int winvmod, long linvmod)
{
    unsigned long *R = flint_malloc((2 * wB) * (sizeof(unsigned long)));
    gf2x_sqr(R, B, wB, lB);
    *lA = 2 * lB - 1;
    *wA = L2W(*lA);

    gf2x_reduce(A, R, *wA, *lA, mod, wmod, lmod, invmod, winvmod, linvmod);
    *lA = FLINT_MIN(lmod - 1, *lA);
    *wA = L2W(*lA);
    gf2x_normalise(A, wA, lA);

    flint_free(R);
}

void padic_poly_get_gf2x(unsigned long **rop, unsigned int *w, long *l,
                         const padic_poly_t op)
{
    if (op->val > 0)
    {
        *l = 0;
        *w = 0;
    }
    else
    {
        long i, j;
        unsigned long bit;

        *l = padic_poly_degree(op) + 1;
        *w = L2W(*l);
        *rop = (unsigned long *) flint_realloc(*rop, *w * sizeof(unsigned long));

        j = -1;
        bit = 1UL;
        for (i = 0; i < *l; i++)
        {
            if (i % (sizeof(unsigned long) * CHAR_BIT) == 0)
            {
                j++;
                bit = 1UL;
                (*rop)[j] = 0;
            }
            if (fmpz_tstbit(op->coeffs + i, 0))
                (*rop)[j] |= bit;
            bit <<= 1;
        }
    }
}

void padic_poly_set_gf2x(padic_poly_t rop, unsigned long *op, unsigned int w, long l)
{
    unsigned int i;
    long k;
    unsigned long bit;

    padic_poly_fit_length(rop, l);
    _padic_poly_set_length(rop, l);
    rop->val = 0;

    i = 0;
    k = 0;
    while (k < l)
    {
        bit = 1UL;
        while (bit != 0 && k < l)
        {
            if (op[i] & bit)
                fmpz_one(rop->coeffs + k);
            else
                fmpz_zero(rop->coeffs + k);
            bit <<= 1;
            k++;
        }
        i++;
    }

    _padic_poly_normalise(rop);
}

/* This really takes place in a finite field of even characteristic */
void _qadic_dense_proot_basis_init_gf2x(qadic_dense_struct *roots, const qadic_dense_ctx_t ctx)
{
    unsigned long *mod = NULL, *invmod = NULL, *r;
    unsigned int wmod, winvmod, w;
    long lmod, linvmod, l;

    padic_poly_get_gf2x(&mod, &wmod, &lmod, ctx->mod);
    padic_poly_get_gf2x(&invmod, &winvmod, &linvmod, ctx->invmod);

    {
        qadic_dense_one(roots + 0, ctx);
    }
    {
        long j = lmod - 1;

        r = (unsigned long *) flint_malloc((2 * wmod) * sizeof(unsigned long));
        r[0] = 2UL;
        w = 1U;
        l = 2L;

        while (--j)
        {
            _qadic_dense_sqr_gf2x(r, &w, &l, r, w, l,
                                  mod, wmod, lmod, invmod, winvmod, linvmod);
        }

        padic_poly_set_gf2x(roots + 1, r, w, l);

        gf2x_clear(mod);
        gf2x_clear(invmod);
        gf2x_clear(r);
    }
}

/* This really takes place in a finite field of even characteristic */
void _qadic_dense_proot_basis_init_char_2(qadic_dense_struct *roots, const qadic_dense_ctx_t ctx)
{
    const fmpz *p = (&ctx->pctx)->p;
    const long d = qadic_dense_ctx_degree(ctx);
    long j;
    qadic_dense_t t, r;
    fmpz_t pow;

    fmpz_init(pow);
    fmpz_pow_ui(pow, p, d - 1);

    padic_poly_init2(t, d);
    _padic_poly_set_length(t, 2);
    fmpz_set_ui(t->coeffs, 0);
    fmpz_set_ui(t->coeffs + 1, 1);
    qadic_dense_pow(t, t, pow, ctx);
    /*j = d;
    while (--j)
    {
        qadic_dense_sqr(t, t, ctx);
        }*/

    qadic_dense_init(r);
    qadic_dense_one(r, ctx);

    j = 0;
    {
        qadic_dense_set(roots, r);
    }
    for (j++; j < *p; j++)
    {
        qadic_dense_mul(r, r, t, ctx);
        qadic_dense_set(roots + j, r);
    }

    qadic_dense_clear(t);
    qadic_dense_clear(r);
}

/* This really takes place in a finite field of small characteristic */
void _qadic_dense_proot_basis_init(qadic_dense_struct *roots, const qadic_dense_ctx_t ctx)
{
    const fmpz *p = (&ctx->pctx)->p;
    const long d = qadic_dense_ctx_degree(ctx);
    long j;
    qadic_dense_t t, r;
    fmpz_t pow;

    fmpz_init(pow);
    fmpz_pow_ui(pow, p, d - 1);

    padic_poly_init2(t, d);
    _padic_poly_set_length(t, 2);
    fmpz_set_ui(t->coeffs, 0);
    fmpz_set_ui(t->coeffs + 1, 1);
    qadic_dense_pow(t, t, pow, ctx);

    qadic_dense_init(r);
    qadic_dense_one(r, ctx);

    j = 0;
    {
        qadic_dense_set(roots, r);
    }
    for (j++; j < *p; j++)
    {
        qadic_dense_mul(r, r, t, ctx);
        qadic_dense_set(roots + j, r);
    }

    qadic_dense_clear(t);
    qadic_dense_clear(r);
}

void qadic_dense_proot_basis_init(qadic_dense_struct **roots, const qadic_dense_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;
    const long p = *(&ctx->pctx)->p;

    if (N != 1)
    {
        printf("Error (qadic_dense_init_proot_basis).  precision different from one\n");
        abort();
    }
    else if (COEFF_IS_MPZ(p))
    {
        printf("Error (qadic_dense_init_proot_basis).  characteristic does not fit in a long\n");
        abort();
    }
    else if (p == 2L)
    {
        long i;
        *roots = flint_malloc(p * sizeof(qadic_dense_struct));
        for (i = 0; i < p; i++)
            qadic_dense_init(*roots + i);
        _qadic_dense_proot_basis_init_gf2x(*roots, ctx);
        /*_qadic_dense_proot_basis_init_char_2(*roots, ctx);*/
    }
    else
    {
        long i;
        *roots = flint_malloc(p * sizeof(qadic_dense_struct));
        for (i = 0; i < p; i++)
            qadic_dense_init(*roots + i);
        _qadic_dense_proot_basis_init(*roots, ctx);
    }
}

void _qadic_dense_proot_gf2x(qadic_dense_t rop, const qadic_dense_t op, const qadic_dense_struct *roots,
                  const qadic_dense_ctx_t ctx)
{
    long k;
    unsigned long *mod = NULL, *invmod = NULL, *sqrt = NULL, *r, *s;
    unsigned int wmod, winvmod, wsqrt, wr, ws;
    long lmod, linvmod, lsqrt, lr, ls;
    unsigned long bit;

    padic_poly_get_gf2x(&mod, &wmod, &lmod, ctx->mod);
    padic_poly_get_gf2x(&invmod, &winvmod, &linvmod, ctx->invmod);
    padic_poly_get_gf2x(&sqrt, &wsqrt, &lsqrt, roots + 1);

    r =  (unsigned long *) flint_malloc((2 * wmod) * sizeof(unsigned long));
    s =  (unsigned long *) flint_malloc((2 * wmod) * sizeof(unsigned long));

    {
        bit = 1UL;
        wr = 0;
        r[wr] = 0;
        for (k = 0; k < op->length; k += 2)
        {
            if (bit == 0)
            {
                bit = 1UL;
                wr++;
                r[wr] = 0;
            }

            if (fmpz_tstbit(op->coeffs + k, 0))
                r[wr] |= bit;

            bit <<= 1;
        }
        lr = (op->length + 1) >> 1;
        wr++;
        gf2x_normalise(r, &wr, &lr);
    }
    {
        bit = 1UL;
        ws = 0;
        s[ws] = 0;
        for (k = 1; k < op->length; k += 2)
        {
            if (bit == 0)
            {
                bit = 1UL;
                ws++;
                s[ws] = 0;
            }

            if (fmpz_tstbit(op->coeffs + k, 0))
                s[ws] |= bit;

            bit <<= 1;
        }
        ls = op->length >> 1;
        ws++;
        gf2x_normalise(s, &ws, &ls);

        gf2x_mul(s, s, ws, sqrt, wsqrt);
        ls = ls + lsqrt - 1;
        ws = L2W(ls);

        gf2x_add(s, s, ws, r, wr);

        gf2x_reduce(r, s, ws, ls, mod, wmod, lmod, invmod, winvmod, linvmod);
        lr = FLINT_MIN(lmod - 1, ls);
        wr = L2W(lr);
        gf2x_normalise(r, &wr, &lr);
    }

    padic_poly_set_gf2x(rop, r, wr, lr);

    flint_free(r);
    flint_free(s);
}

void _qadic_dense_proot(qadic_dense_t rop, const qadic_dense_t op, const qadic_dense_struct *roots,
                  const qadic_dense_ctx_t ctx)
{
    const long d = qadic_dense_ctx_degree(ctx);
    const long p = *(&ctx->pctx)->p;
    long j, k;
    qadic_dense_t t;

    qadic_dense_init(t);
    padic_poly_fit_length(t, d);

    j = 0;
    {
        qadic_dense_zero(rop);
        for (k = 0; p*k + j < op->length; k++)
            fmpz_set(rop->coeffs + k, op->coeffs + p*k + j);
        _padic_poly_set_length(rop, k);
        rop->val = op->val;
        padic_poly_canonicalise(rop, &p);
        _padic_poly_normalise(rop);
    }
    for (j = 1; j < p; j++)
    {
        for (k = 0; p*k + j < op->length; k++)
            fmpz_set(t->coeffs + k, op->coeffs + p*k + j);
        _padic_poly_set_length(t, k);
        t->val = op->val;
        padic_poly_canonicalise(t, &p);
        _padic_poly_normalise(t);

        qadic_dense_mul(t, t, roots + j, ctx);

        qadic_dense_add(rop, rop, t, ctx);
    }

    qadic_dense_clear(t);
}

void qadic_dense_proot(qadic_dense_t rop, const qadic_dense_t op,
                 const qadic_dense_struct *roots,
                 const qadic_dense_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;
    const long p = *(&ctx->pctx)->p;

    if (N != 1)
    {
        printf("ERROR (qadic_dense_proot).  not implemented.\n");
        abort();
    }
    else if (COEFF_IS_MPZ(p))
    {
        printf("Error (qadic_dense_proot).  characteristic does not fit in a long\n");
        abort();
    }
    else if (qadic_dense_is_zero(op))
    {
        qadic_dense_zero(rop);
        return;
    }
#if 1
    else if (p == 2L)
    {
        if (rop == op)
        {
            qadic_dense_t t;
            qadic_dense_init(t);
            qadic_dense_set(t, op);
            _qadic_dense_proot_gf2x(rop, t, roots, ctx);
            qadic_dense_clear(t);
        }
        else
        {
            _qadic_dense_proot_gf2x(rop, op, roots, ctx);
        }
    }
#endif
    else
    {
        if (rop == op)
        {
            qadic_dense_t t;
            qadic_dense_init(t);
            qadic_dense_set(t, op);
            _qadic_dense_proot(rop, t, roots, ctx);
            qadic_dense_clear(t);
        }
        else
        {
            _qadic_dense_proot(rop, op, roots, ctx);
        }
    }
}

void _padic_poly_add_char_2(fmpz *rop, long *val,
                          const fmpz *op1, long val1, long len1,
                          const fmpz *op2, long val2, long len2,
                          const padic_ctx_t ctx)
{
    const long len = FLINT_MAX(len1, len2);

    *val = FLINT_MIN(val1, val2);

    if (val1 == val2)
    {
        _fmpz_poly_add(rop, op1, len1, op2, len2);
    }
    else  /* => (op1 != op2) */
    {
        long exp;

        if (val1 < val2)  /* F := p^g (G + p^{h-g} H) */
        {
            exp = val2 - val1;

            if (rop == op1)
            {
                _fmpz_vec_zero(rop + len1, len2 - len1);
                _fmpz_vec_scalar_addmul_si_2exp(rop, op2, len2, 1L, exp);
            }
            else
            {
                _fmpz_vec_scalar_mul_2exp(rop, op2, len2, exp);
                _fmpz_poly_add(rop, op1, len1, rop, len2);
            }
        }
        else  /* F := p^h (p^{g-h} G + H) */
        {
            exp = val1 - val2;

            if (rop == op2)
            {
                _fmpz_vec_zero(rop + len2, len1 - len2);
                _fmpz_vec_scalar_addmul_si_2exp(rop, op1, len1, 1L, exp);
            }
            else
            {
                _fmpz_vec_scalar_mul_2exp(rop, op1, len1, exp);
                _fmpz_poly_add(rop, rop, len1, op2, len2);
            }
        }
    }

    _padic_poly_canonicalise(rop, val, len, ctx->p);
}

void _padic_poly_add_no_mod(fmpz *rop, long *val, 
                     const fmpz *op1, long val1, long len1, 
                     const fmpz *op2, long val2, long len2, 
                     const padic_ctx_t ctx)
{
    const long len = FLINT_MAX(len1, len2);
    fmpz_t pow;
    int alloc;

    *val = FLINT_MIN(val1, val2);

    alloc = _padic_ctx_pow_ui(pow, ctx->N - *val, ctx);

    if (val1 == val2)
    {
        _fmpz_poly_add(rop, op1, len1, op2, len2);
    }
    else  /* => (op1 != op2) */
    {
        fmpz_t x;

        fmpz_init(x);
        if (val1 < val2)  /* F := p^g (G + p^{h-g} H) */
        {
            fmpz_pow_ui(x, ctx->p, val2 - val1);

            if (rop == op1)
            {
                _fmpz_vec_zero(rop + len1, len2 - len1);
                _fmpz_vec_scalar_addmul_fmpz(rop, op2, len2, x);
            }
            else
            {
                _fmpz_vec_scalar_mul_fmpz(rop, op2, len2, x);
                _fmpz_poly_add(rop, op1, len1, rop, len2);
            }
        }
        else  /* F := p^h (p^{g-h} G + H) */
        {
            fmpz_pow_ui(x, ctx->p, val1 - val2);

            if (rop == op2)
            {
                _fmpz_vec_zero(rop + len2, len1 - len2);
                _fmpz_vec_scalar_addmul_fmpz(rop, op1, len1, x);
            }
            else
            {
                _fmpz_vec_scalar_mul_fmpz(rop, op1, len1, x);
                _fmpz_poly_add(rop, rop, len1, op2, len2);
            }
        }
        fmpz_clear(x);
    }

    if (alloc)
        fmpz_clear(pow);

    _padic_poly_canonicalise(rop, val, len, ctx->p);
}

void padic_poly_add_no_mod(padic_poly_t f, 
                    const padic_poly_t g, const padic_poly_t h, 
                    const padic_ctx_t ctx)
{
    const long lenG = g->length;
    const long lenH = h->length;
    const long lenF = FLINT_MAX(lenG, lenH);

    if (lenG == 0 && lenH == 0)
    {
        padic_poly_zero(f);
        return;
    }

    padic_poly_fit_length(f, lenF);

    _padic_poly_add_no_mod(f->coeffs, &(f->val), g->coeffs, g->val, lenG, 
                                          h->coeffs, h->val, lenH, ctx);

    _padic_poly_set_length(f, lenF);
    _padic_poly_normalise(f);
}

void _qadic_dense_frobenius_teichmuller_char_2(qadic_dense_t rop, const qadic_dense_t op,
                                 const qadic_dense_ctx_t ctx)
{
    const long len = padic_poly_length(op);
    long i;


    if (len == 0)
    {
        qadic_dense_zero(rop);
    }
    else
    {
        padic_poly_fit_length(rop, 2*len - 1);
        _fmpz_vec_zero(rop->coeffs, 2*len - 1);
        i = 0;
        {
            fmpz_set(rop->coeffs, op->coeffs);
        }
        for (i++; i < len; i++)
        {
            fmpz_zero(rop->coeffs + 2*i - 1);
            fmpz_set(rop->coeffs + 2*i, op->coeffs + i);
        }
        _padic_poly_set_length(rop, 2*len - 1);
        rop->val = op->val;

        qadic_dense_reduce(rop, ctx);
    }
}

void _qadic_dense_frobenius_teichmuller(qadic_dense_t rop, const qadic_dense_t op,
                                 const qadic_dense_ctx_t ctx)
{
    const long len = padic_poly_length(op);
    const long p = *(&ctx->pctx)->p;
    long i;

    if (len == 0)
    {
        qadic_dense_zero(rop);
    }
    else
    {
        padic_poly_fit_length(rop, p*(len - 1) + 1);
        _fmpz_vec_zero(rop->coeffs, p*(len - 1) + 1);
        for (i = 0; i < len; i++)
        {
            fmpz_set(rop->coeffs + p*i, op->coeffs + i);
        }
        _padic_poly_set_length(rop, p*(len - 1) + 1);
        rop->val = op->val;

        qadic_dense_reduce(rop, ctx);
    }
}

void qadic_dense_frobenius_teichmuller(qadic_dense_t rop, const qadic_dense_t op,
                                 const qadic_dense_ctx_t ctx)
{
    if (rop == op)
    {
        qadic_dense_t t;
        qadic_dense_set(t, op);
        _qadic_dense_frobenius_teichmuller(rop, t, ctx);
        qadic_dense_clear(t);
    }
    else
    {
        _qadic_dense_frobenius_teichmuller(rop, op, ctx);
    }
}

void _qadic_dense_artin_schreier_root_ii(qadic_dense_t x, const qadic_dense_t alpha,
                                   const qadic_dense_t beta, const qadic_dense_t gamma,
                                   const qadic_dense_struct *roots,
                                   long m, const long *b, const qadic_dense_ctx_struct *qctx)
{
    if (m == 1)
    {
        /* This occurs in F_2[X] */
#if 0
        qadic_dense_t alpha2, gamma2;

        qadic_dense_init(alpha2);
        qadic_dense_init(gamma2);

        qadic_dense_set(alpha2, alpha);
        qadic_dense_set(gamma2, gamma);

        padic_poly_reduce(alpha2, &qctx->pctx);
        padic_poly_reduce(gamma2, &qctx->pctx);

        qadic_dense_inv(x, alpha2, qctx);
        qadic_dense_mul(x, x, gamma2, qctx);
        qadic_dense_neg(x, x, qctx);
        qadic_dense_proot(x, x, roots, qctx);

        qadic_dense_clear(alpha2);
        qadic_dense_clear(gamma2);
#else
        unsigned long *mod = NULL, *invmod = NULL, *sqrt = NULL, *alpha2 = NULL, *gamma2 = NULL, *r, *s;
        unsigned int wmod, winvmod, wsqrt, wa2, wc2, wr, ws;
        long lmod, linvmod, lsqrt, la2, lc2, lr, ls;

        padic_poly_get_gf2x(&mod, &wmod, &lmod, qctx->mod);
        padic_poly_get_gf2x(&invmod, &winvmod, &linvmod, qctx->invmod);
        padic_poly_get_gf2x(&sqrt, &wsqrt, &lsqrt, roots + 1);

        padic_poly_get_gf2x(&alpha2, &wa2, &la2, alpha);
        padic_poly_get_gf2x(&gamma2, &wc2, &lc2, gamma);

        r = flint_malloc((2 * wmod) * sizeof(unsigned long));
        s = flint_malloc((2 * wmod) * sizeof(unsigned long));

        gf2x_invmod(r, alpha2, wa2, la2, mod, wmod, lmod);
        lr = lmod - 1;
        wr = L2W(lr);
        gf2x_normalise(r, &wr, &lr);

        gf2x_mul(r, r, wr, gamma2, wc2);
        lr = lr + lc2 - 1;
        wr = L2W(lr);
        gf2x_reduce(s, r, wr, lr, mod, wmod, lmod, invmod, winvmod, linvmod);
        ls = FLINT_MIN(lmod - 1, lr);
        ws = L2W(ls);
        gf2x_normalise(s, &ws, &ls);

        _qadic_dense_sqrt_gf2x(s, &ws, &ls, s, ws, ls, sqrt, wsqrt, lsqrt,
                               mod, wmod, lmod, invmod, winvmod, linvmod);
        padic_poly_set_gf2x(x, s, ws, ls);

        flint_free(r);
        flint_free(s);
        gf2x_clear(mod);
        gf2x_clear(invmod);
        gf2x_clear(sqrt);
        gf2x_clear(alpha2);
        gf2x_clear(gamma2);
#endif
    }
    else
    {
        const qadic_dense_ctx_struct *qctxi, *qctxj;
        long *a, h, i, j, n;
        qadic_dense_t alpha2, beta2, gamma2, gamma3, Delta2;

        /* Number of steps */
        n = FLINT_CLOG2(b[0]) + 1;

        /* Precision needed at each step */
        a = flint_malloc(n * sizeof(long));

        i = 0;
        for (a[i = 0] = b[0]; a[i] > 1; i++)
        {
            a[i + 1] = (a[i] + 1) / 2;
        }
        qadic_dense_init(alpha2);
        qadic_dense_init(beta2);
        qadic_dense_init(gamma2);
        qadic_dense_init(gamma3);
        qadic_dense_init(Delta2);

        /* Lifting */
        i = n - 1;
        j = m - 1;
        {
            /* This occurs in F_2[X] */
#if 0
            qctxi = qctx + j;

            qadic_dense_set(alpha2, alpha);
            qadic_dense_set(gamma2, gamma);

            padic_poly_reduce(alpha2, &qctxi->pctx);
            padic_poly_reduce(gamma2, &qctxi->pctx);

            qadic_dense_inv(x, alpha2, qctxi);
            qadic_dense_mul(x, x, gamma2, qctxi);
            qadic_dense_neg(x, x, qctxi);
            qadic_dense_proot(x, x, roots, qctxi);
#else
            unsigned long *mod = NULL, *invmod = NULL, *sqrt = NULL, *a2 = NULL, *c2 = NULL, *r, *s;
            unsigned int wmod, winvmod, wsqrt, wa2, wc2, wr, ws;
            long lmod, linvmod, lsqrt, la2, lc2, lr, ls;

            qctxi = qctx + j;

            padic_poly_get_gf2x(&mod, &wmod, &lmod, qctxi->mod);
            padic_poly_get_gf2x(&invmod, &winvmod, &linvmod, qctxi->invmod);
            padic_poly_get_gf2x(&sqrt, &wsqrt, &lsqrt, roots + 1);

            padic_poly_get_gf2x(&a2, &wa2, &la2, alpha);
            padic_poly_get_gf2x(&c2, &wc2, &lc2, gamma);

            gf2x_normalise(a2, &wa2, &la2);
            gf2x_normalise(c2, &wc2, &lc2);

            r = flint_malloc((2 * wmod) * sizeof(unsigned long));
            s = flint_malloc((2 * wmod) * sizeof(unsigned long));

            gf2x_invmod(r, a2, wa2, la2, mod, wmod, lmod);
            lr = lmod - 1;
            wr = L2W(lr);
            gf2x_normalise(r, &wr, &lr);

            gf2x_mul(r, r, wr, c2, wc2);
            lr = lr + lc2 - 1;
            wr = L2W(lr);
            gf2x_reduce(s, r, wr, lr, mod, wmod, lmod, invmod, winvmod, linvmod);
            ls = FLINT_MIN(lmod - 1, lr);
            ws = L2W(ls);
            gf2x_normalise(s, &ws, &ls);

            _qadic_dense_sqrt_gf2x(s, &ws, &ls, s, ws, ls, sqrt, wsqrt, lsqrt,
                                   mod, wmod, lmod, invmod, winvmod, linvmod);
            padic_poly_set_gf2x(x, s, ws, ls);

            flint_free(r);
            flint_free(s);
            gf2x_clear(mod);
            gf2x_clear(invmod);
            gf2x_clear(sqrt);
            gf2x_clear(a2);
            gf2x_clear(c2);
#endif
        }
        for (i--; i >= 0; i--)
        {
            h = j;
            j--;
            qctxi--;
            if (b[j] != a[i])
            {
                j--;
                qctxi--;
            }
            if (b[j + 1] == a[i] - a[i + 1])
                qctxj = qctx + j + 1;
            else
                qctxj = qctx + j + 2;

            qadic_dense_set(alpha2, alpha);
            qadic_dense_set(beta2, beta);
            qadic_dense_set(gamma2, gamma);

            padic_poly_reduce(alpha2, &qctxi->pctx);
            padic_poly_reduce(beta2, &qctxi->pctx);
            padic_poly_reduce(gamma2, &qctxi->pctx);

            qadic_dense_frobenius_teichmuller(Delta2, x, qctxi);
            padic_poly_mul(gamma3, Delta2, alpha2, &qctxi->pctx);
            padic_poly_mul(Delta2, beta2, x, &qctxi->pctx);
            padic_poly_add(gamma3, gamma3, Delta2, &qctxi->pctx);
            qadic_dense_reduce(gamma3, qctxi);
            qadic_dense_add(gamma3, gamma3, gamma2, qctxi);
            if (!qadic_dense_is_zero(gamma3))
            {
                gamma3->val -= a[i + 1];
            }

            padic_poly_reduce(alpha2, &qctxj->pctx);
            padic_poly_reduce(beta2, &qctxj->pctx);
            padic_poly_reduce(gamma3, &qctxj->pctx);

            if (b[h] == a[i] - a[i + 1])
                _qadic_dense_artin_schreier_root_ii(Delta2, alpha2, beta2, gamma3, roots, m - h, b + h, qctx + h);
            else
                _qadic_dense_artin_schreier_root_ii(Delta2, alpha2, beta2, gamma3, roots, m - h - 1, b + h + 1, qctx + h + 1);
            Delta2->val += a[i + 1];

            qadic_dense_add(x, x, Delta2, qctxi);
        }

        qadic_dense_clear(alpha2);
        qadic_dense_clear(beta2);
        qadic_dense_clear(gamma2);
        qadic_dense_clear(gamma3);
        qadic_dense_clear(Delta2);

        flint_free(a);
    }
}

void qadic_dense_artin_schreier_root_ii(qadic_dense_t x, const qadic_dense_t alpha,
                                   const qadic_dense_t beta, const qadic_dense_t gamma,
                                   const qadic_dense_struct *roots,
                                   const qadic_dense_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (N == 1)
    {
        qadic_dense_ctx_t qctx;
        qadic_dense_t alpha2, gamma2;

        qadic_dense_ctx_init_reduce(qctx, ctx, 1);

        qadic_dense_init(alpha2);
        qadic_dense_init(gamma2);

        qadic_dense_set(alpha2, alpha);
        qadic_dense_set(gamma2, gamma);

        padic_poly_reduce(alpha2, &qctx->pctx);
        padic_poly_reduce(gamma2, &qctx->pctx);

        qadic_dense_inv(x, alpha2, qctx);
        qadic_dense_mul(x, x, gamma2, qctx);
        qadic_dense_neg(x, x, qctx);
        qadic_dense_proot(x, x, roots, qctx);

        qadic_dense_clear(alpha2);
        qadic_dense_clear(gamma2);

        qadic_dense_ctx_clear(qctx);
    }
    else
    {
        long *b, i, n, m;
        qadic_dense_ctx_struct *qctx;

        /* Number of steps */
        n = FLINT_CLOG2(N) + 1;

        /* Compute precisions needed later */
        b = flint_malloc((2 * n) * sizeof(long));

        i = 0;
        {
            b[i] = N;
        }
        for (; !(b[i] & 1L); i++)
        {
            b[i + 1] = b[i] >> 1;
        }
        for (; b[i] > 2; i += 2)
        {
            b[i + 1] = (b[i] >> 1) + 1;
            b[i + 2] = b[i] >> 1;
        }
        if (b[i] == 2L)
        {
            b[i + 1] = 1L;
            i++;
        }

        /* Contexts */
        m = i + 1;
        qctx = flint_malloc(m * sizeof(qadic_dense_ctx_struct));

        i = 0;
        {
            qadic_dense_ctx_init_reduce(qctx + i, ctx, b[i]);
        }
        for (i++; i < m; i++)
        {
            qadic_dense_ctx_init_reduce(qctx + i, qctx + i - 1, b[i]);
        }

        _qadic_dense_artin_schreier_root_ii(x, alpha, beta, gamma, roots, m, b, qctx);

        for (i = 0; i < m; i++)
            qadic_dense_ctx_clear(qctx + i);

        flint_free(qctx);

        flint_free(b);
    }
}

/* Newton lift for the modified modular polynomial \Gamma_2 */
/* \Gamma_2 = (X + 2 Y + 8 X Y)^2 + Y + 4 X Y */
/* k = 0, M = N' */
void _qadic_dense_gen_newton_lift_ii_gamma_2(qadic_dense_t rop, const qadic_dense_t op,
                                      const qadic_dense_struct *roots,
                                      long m, const long *b, const qadic_dense_ctx_struct *qctx)
{
    if (m == 1)
    {
        qadic_dense_set(rop, op);
    }
    else if (*(&qctx->pctx)->p == 2L)
    {
        const qadic_dense_ctx_struct *qctxi, *qctxk;
        long *a, i, j, k, n;
        qadic_dense_t Delta_x, Delta_y, Delta, V, y, xy;

        /* Number of steps */
        n = FLINT_CLOG2(b[0]) + 1;

        /* Precision needed at each step */
        a = flint_malloc(n * sizeof(long));

        i = 0;
        for (a[i = 0] = b[0]; a[i] > 1; i++)
        {
            a[i + 1] = (a[i] + 1) / 2;
        }

        qadic_dense_init(Delta_x);
        qadic_dense_init(Delta_y);
        qadic_dense_init(Delta);
        qadic_dense_init(V);
        qadic_dense_init(y);
        qadic_dense_init(xy);

        /* Lifting */
        i = n - 1;
        j = m - 1;
        {
            qctxi = qctx + j;

            qadic_dense_set(rop, op);
        }
        for (i--; i >= 0; i--)
        {
            /* Previous step precision */
            /* This is the precision needed for the Artin--Schreier subroutine as well */
            k = j;
            qctxk = qctxi;
            /* Current step precision */
            j--;
            qctxi--;
            if (b[j] != a[i])
            {
                j--;
                qctxi--;
            }

            /* y */
            qadic_dense_frobenius_teichmuller(y, rop, qctxi);
            /* 2 y */
            qadic_dense_set(Delta_x, y);
            Delta_x->val += 1;
            padic_poly_reduce(Delta_x, &qctxi->pctx);
            /* 2 x y */
            qadic_dense_mul(xy, rop, Delta_x, qctxi);
            /* x + 2 y */
            qadic_dense_add(Delta, rop, Delta_x, qctxi);
            /* 8 x y */
            qadic_dense_set(Delta_x, xy);
            Delta_x->val += 2;
            padic_poly_reduce(Delta_x, &qctxi->pctx);
            /* x + 2 y + 8 x y */
            qadic_dense_add(Delta, Delta, Delta_x, qctxi);

            /* V */
            /* (x + 2 y + 8 x y)² */
            qadic_dense_sqr(V, Delta, qctxi);
            /* (x + 2 y + 8 x y)² + y */
            padic_poly_add_no_mod(V, V, y, &qctxi->pctx);
            /* 4 x y */
            xy->val += 1;
            /* (x + 2 y + 8 x y)² + y + 4 x y */
            padic_poly_add_no_mod(V, V, xy, &qctxi->pctx);
            padic_poly_reduce(V, &qctxi->pctx);
            /* (x + 2 y + 8 x y)² + y + 4 x y / 2^N' */
            V->val -= a[i + 1];

            /* Delta_x */
            /* 2 (x + 2 y + 8 x y) */
            Delta->val += 1;
            padic_poly_reduce(Delta, &qctxk->pctx);
            /* 1 */
            qadic_dense_one(xy, qctxk);
            /* 8 y */
            qadic_dense_set(Delta_x, y);
            Delta_x->val += 3;
            padic_poly_reduce(Delta_x, &qctxk->pctx);
            /* 1 + 8 y */
            qadic_dense_add(xy, xy, Delta_x, qctxk);
            /* 2 (x + 2 y + 8 x y) (1 + 8 y) */
            qadic_dense_mul(Delta_x, Delta, xy, qctxk);
            /* 4 y */
            y->val += 2;
            padic_poly_reduce(y, &qctxk->pctx);
            /* 2 (x + 2 y + 8 x y) (1 + 8 y) + 4 y*/
            qadic_dense_add(Delta_x, Delta_x, y, qctxk);

#if 0
            /* Delta_y */
            /* 4 (x + 2 y + 8 x y) */
            Delta->val += 1;
            padic_poly_reduce(Delta, &qctxk->pctx);
            /* 1 + 4 (x + 2 y + 8 x y) */
            qadic_dense_one(xy, qctxk);
            qadic_dense_add(Delta_y, Delta, xy, qctxk);
            /* 4 x */
            qadic_dense_set(y, rop);
            y->val += 2;
            padic_poly_reduce(y, &qctxk->pctx);
            /* 1 + 4 x */
            qadic_dense_add(xy, xy, y, qctxk);
            /* (1 + 4 (x + 2 y + 8 x y)) (1 + 4 x) */
            qadic_dense_mul(Delta_y, Delta_y, xy, qctxk);
#else
            /* Delta_y */
            /* 4 (x + 2 y + 8 x y) */
            Delta->val += 1;
            padic_poly_reduce(Delta, &qctxk->pctx);
            /* 4 x */
            qadic_dense_set(y, rop);
            y->val += 2;
            padic_poly_reduce(y, &qctxk->pctx);
            /* 16 x * (x + 2 y + 8 x y) */
            qadic_dense_mul(Delta_y, Delta, y, qctxk);
            /* 4 (x + 2 y + 8 x y) + 16 x * (x + 2 y + 8 x y) */
            padic_poly_add_no_mod(Delta_y, Delta_y, Delta, &qctxk->pctx);
            /* 4 x + 4 (x + 2 y + 8 x y) + 16 x * (x + 2 y + 8 x y) */
            padic_poly_add_no_mod(Delta_y, Delta_y, y, &qctxk->pctx);
            /* 1 + 4 x + 4 (x + 2 y + 8 x y) + 16 x * (x + 2 y + 8 x y) */
            qadic_dense_one(xy, qctxk);
            padic_poly_add_no_mod(Delta_y, Delta_y, xy, &qctxk->pctx);
            padic_poly_reduce(Delta_y, &qctxk->pctx);
#endif
            /* Delta */
            _qadic_dense_artin_schreier_root_ii(Delta, Delta_y, Delta_x, V, roots, m - k, b + k, qctx + k);
            Delta->val += a[i + 1];
            /*padic_poly_reduce(Delta, &qctxi->pctx);*/

            /* rop */
            qadic_dense_add(rop, rop, Delta, qctxi);

#if DEBUG >= 1
            /* y */
            qadic_dense_frobenius_teichmuller(y, rop, qctxi);
            /* 2 y */
            qadic_dense_set(Delta_x, y);
            Delta_x->val += 1;
            padic_poly_reduce(Delta_x, &qctxi->pctx);
            /* 2 x y */
            qadic_dense_mul(xy, rop, Delta_x, qctxi);
            /* x + 2 y */
            qadic_dense_add(Delta, rop, Delta_x, qctxi);
            /* 8 x y */
            qadic_dense_set(Delta_x, xy);
            Delta_x->val += 2;
            padic_poly_reduce(Delta_x, &qctxi->pctx);
            /* x + 2 y + 8 x y */
            qadic_dense_add(Delta, Delta, Delta_x, qctxi);

            /* V */
            /* (x + 2 y + 8 x y)² */
            qadic_dense_sqr(V, Delta, qctxi);
            /* (x + 2 y + 8 x y)² + y */
            qadic_dense_add(V, V, y, qctxi);
            /* 4 x y */
            xy->val += 1;
            padic_poly_reduce(xy, &qctxi->pctx);
            /* (x + 2 y + 8 x y)² + y + 4 x y */
            qadic_dense_add(V, V, xy, qctxi);
            printf("Gamma [%ld] ", a[i]), _fmpz_vec_print(V->coeffs, V->length), printf("\n");
#endif
        }

        qadic_dense_clear(Delta_x);
        qadic_dense_clear(Delta_y);
        qadic_dense_clear(Delta);
        qadic_dense_clear(V);
        qadic_dense_clear(y);
        qadic_dense_clear(xy);
    }
    else
    {
        printf("ERROR (_padic_poly_teichmuller).  odd characteristic not implemented.\n");
        abort();
    }
}

/* Newton lift for the modified modular polynomial \Gamma_2 */
/* \Gamma_2 = (X + 2 Y + 8 X Y)^2 + Y + 4 X Y */
/* k = 0, M = N' */
void qadic_dense_gen_newton_lift_ii_gamma_2(qadic_dense_t rop, const qadic_dense_t op,
                                      const qadic_dense_struct *roots,
                                      const qadic_dense_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (N == 1)
    {
        qadic_dense_set(rop, op);
    }
    else if (*(&ctx->pctx)->p == 2L)
    {
        long *b, i, m, n;
        qadic_dense_ctx_struct *qctx;

        /* Number of steps */
        n = FLINT_CLOG2(N) + 1;

        /* Compute precisions needed later */
        b = flint_malloc((2 * n) * sizeof(long));

        i = 0;
        /* Precision to reach */
        {
            b[i] = N;
        }
        /* As long as the previous precision is even,
           only one additional precision will be needed for the next steps */
        for (; !(b[i] & 1L); i++)
        {
            b[i + 1] = b[i] >> 1;
        }
        /* Once a precision is odd, floor and ceil of its half will be needed */
        for (; b[i] > 2; i += 2)
        {
            b[i + 1] = (b[i] >> 1) + 1;
            b[i + 2] = b[i] >> 1;
        }
        /* 2 needs to be treated separately to avoid repetition */
        if (b[i] == 2L)
        {
            b[i + 1] = 1L;
            i++;
        }

        /* Contexts */
        m = i + 1;
        qctx = flint_malloc(m * sizeof(qadic_dense_ctx_struct));

        i = 0;
        {
            qadic_dense_ctx_init_reduce(qctx + i, ctx, b[i]);
        }
        for (i++; i < m; i++)
        {
            qadic_dense_ctx_init_reduce(qctx + i, qctx + i - 1, b[i]);
        }

        _qadic_dense_gen_newton_lift_ii_gamma_2(rop, op, roots, m, b, qctx);

        for (i = 0; i < n; i++)
            qadic_dense_ctx_clear(qctx + i);

        flint_free(qctx);

        flint_free(b);
    }
    else
    {
        printf("ERROR (_padic_poly_teichmuller).  odd characteristic not implemented.\n");
        abort();
    }
}

/* Newton lift for the modified modular polynomial \Lambda_2 */
/* \Lambda_2 = Y^2 (1 + X)^2 - 4 X */
/* k = 3 */
void _qadic_dense_gen_newton_lift_ii_lambda_2(qadic_dense_t rop, const qadic_dense_t op,
                                      const qadic_dense_struct *roots,
                                      long N, const qadic_dense_ctx_t ctx)
{
    if (N < 2 * 3 + 1 + 1)
    {
        qadic_dense_set(rop, op);
    }
    else if (*(&ctx->pctx)->p == 2L)
    {
        qadic_dense_ctx_struct *qctxi, *qctxj, *qctxk;

        long *a, i, n;
        qadic_dense_t Delta_x, Delta_y, Delta, V, xy;

        long M, *b, l, m;
        qadic_dense_ctx_struct *qctx;

        /* Number of steps */
        n = FLINT_CLOG2(N - 2 * 3) + 1;

        /* Precision needed at each step */
        a = flint_malloc(n * sizeof(long));

        i = 0;
        for (a[i = 0] = N; a[i] > 2 * 3 + 1; i++)
        {
            a[i + 1] = (a[i] + 1) / 2 + 3;
        }
        /* i = n - 1 */

        qadic_dense_init(Delta_x);
        qadic_dense_init(Delta_y);
        qadic_dense_init(Delta);
        qadic_dense_init(V);
        qadic_dense_init(xy);

        qctxi = flint_malloc(sizeof(qadic_dense_ctx_struct));
        qctxk = flint_malloc(sizeof(qadic_dense_ctx_struct));

        n = (FLINT_CLOG2(N) + 1);
        b = flint_malloc((2 * n) * sizeof(long));
        qctx = flint_malloc((2 * n) * sizeof(qadic_dense_ctx_struct));

        /* Lifting */
        {
            qadic_dense_ctx_init_reduce(qctxi, ctx, a[i]);
            qadic_dense_ctx_init_reduce(qctxk, ctx, a[i]);

            qadic_dense_set(rop, op);
        }
        for (i--; i >= 0; i--)
        {
            /* Compute precisions for Artin--Schreier subroutine */
            M = (a[i] + 1)/2;

            l = 0;
            {
                b[l] = M;
            }
            for (; !(b[l] & 1L); l++)
            {
                b[l + 1] = b[l] >> 1;
            }
            for (; b[l] > 2; l += 2)
            {
                b[l + 1] = (b[l] >> 1) + 1;
                b[l + 2] = b[l] >> 1;
            }
            if (b[l] == 2L)
            {
                b[l + 1] = 1L;
                l++;
            }

            /* Contexts for Artin--Schreier subroutine */
            m = l + 1;
            l = 0;
            {
                qadic_dense_ctx_init_reduce(qctx + l, ctx, b[l]);
            }
            for (l++; l < m; l++)
            {
                qadic_dense_ctx_init_reduce(qctx + l, qctx + l - 1, b[l]);
            }

            /* Context for previous step */
            qctxj = qctxk;
            qctxk = qctxi;
            qctxi = qctxj;

            /* Context for current step */
            qadic_dense_ctx_clear(qctxi);
            qadic_dense_ctx_init_reduce(qctxi, ctx, a[i]);

            /* y */
            qadic_dense_frobenius_teichmuller(Delta_x, rop, qctxi);
            /* 1 */
            qadic_dense_one(Delta_y, qctxi);
            /* 1 + x */
            qadic_dense_add(Delta_y, Delta_y, rop, qctxi);
            /* y (1 + x) */
            qadic_dense_mul(xy, Delta_x, Delta_y, qctxi);

            /* V */
            /* 4 x */
            qadic_dense_set(Delta, rop);
            Delta->val += 2;
            padic_poly_reduce(Delta, &qctxi->pctx);
            /* y² (1 + x)² */
            qadic_dense_sqr(V, xy, qctxi);
            /* y² (1 + x)² - 4 x */
            qadic_dense_sub(V, V, Delta, qctxi);
            /* y² (1 + x)² - 4 x / 2^N'*/
            V->val -= a[i + 1];

            /* y (1 + x) */
            padic_poly_reduce(xy, &qctxk->pctx);

            /* Delta_y */
            /* 1 + x */
            padic_poly_reduce(Delta_y, &qctxk->pctx);
            /* y (1 + x)² */
            qadic_dense_mul(Delta_y, Delta_y, xy, qctxk);
            /* 2 y (1 + x)² */
            Delta_y->val += 1;
            /* 2 y (1 + x)² / 2^k */
            Delta_y->val -= 3;

            /* Delta_x */
            /* y */
            padic_poly_reduce(Delta_x, &qctxk->pctx);
            /* y² (1 + x) */
            qadic_dense_mul(Delta_x, Delta_x, xy, qctxk);
            /* 2 */
            qadic_dense_set_ui(xy, 2UL, qctxk);
            /* y² (1 + x) - 2 */
            qadic_dense_sub(Delta_x, Delta_x, xy, qctxk);
            /* 2 y² (1 + x) - 4 */
            Delta_x->val += 1;
            /* 2 y² (1 + x) - 4 / 2^k */
            Delta_x->val -= 3;

            /* Delta */
            _qadic_dense_artin_schreier_root_ii(Delta, Delta_y, Delta_x, V, roots, m, b, qctx);
            Delta->val += M;

            /* rop */
            qadic_dense_add(rop, rop, Delta, qctxi);
            padic_poly_reduce(rop, &qctxi->pctx);

#if DEBUG >= 1
            qadic_dense_frobenius_teichmuller(Delta_x, rop, qctxi);
            /* 1 */
            qadic_dense_one(Delta_y, qctxi);
            /* 1 + x */
            qadic_dense_add(Delta_y, Delta_y, rop, qctxi);
            /* y (1 + x) */
            qadic_dense_mul(xy, Delta_x, Delta_y, qctxi);

            /* V */
            /* 4 x */
            qadic_dense_set(Delta, rop);
            Delta->val += 2;
            padic_poly_reduce(Delta, &qctxi->pctx);
            /* y² (1 + x)² */
            qadic_dense_sqr(V, xy, qctxi);
            /* y² (1 + x)² - 4 x */
            qadic_dense_sub(V, V, Delta, qctxi);
            padic_poly_reduce(V, &qctxi->pctx);
            printf("Lambda[%ld] = ", a[i]), qadic_dense_print_pretty(V, qctxi) ,printf("\n");
#endif
            for (l = 0; l < m; l++)
                qadic_dense_ctx_clear(qctx + l);
        }

        qadic_dense_clear(Delta_x);
        qadic_dense_clear(Delta_y);
        qadic_dense_clear(Delta);
        qadic_dense_clear(V);
        qadic_dense_clear(xy);

        qadic_dense_ctx_clear(qctxi);
        qadic_dense_ctx_clear(qctxk);
        flint_free(qctxi);
        flint_free(qctxk);

        flint_free(b);
        flint_free(qctx);
    }
    else
    {
        printf("ERROR (_padic_poly_teichmuller).  odd characteristic not implemented.\n");
        abort();
    }
}

/* Newton lift for the modified modular polynomial \Lambda_2 */
/* \Lambda_2 = Y^2 (1 + X)^2 - 4 X */
/* k = 3 */
void qadic_dense_gen_newton_lift_ii_lambda_2(qadic_dense_t rop, const qadic_dense_t op,
                                      const qadic_dense_struct *roots,
                                      const qadic_dense_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (N < 2 * 3 + 1 + 1)
    {
        qadic_dense_set(rop, op);
    }
    else if (*(&ctx->pctx)->p == 2L)
    {
        _qadic_dense_gen_newton_lift_ii_lambda_2(rop, op, roots, N, ctx);
    }
    else
    {
        printf("ERROR (_padic_poly_teichmuller).  odd characteristic not implemented.\n");
        abort();
    }
}

/* Newton lift for the modified modular polynomial \Lambda'_2 */
/* \Lambda'_2 = 4 X Y^2 - (1 + X)^2 */
/* k = 3 */
void _qadic_dense_gen_newton_lift_ii_invlambda_2(qadic_dense_t rop, const qadic_dense_t op,
                                      const qadic_dense_struct *roots,
                                      long N, const qadic_dense_ctx_t ctx)
{
    if (N < 2 * 3 + 1 + 1)
    {
        qadic_dense_set(rop, op);
    }
    else if (*(&ctx->pctx)->p == 2L)
    {
        qadic_dense_ctx_struct *qctxi, *qctxj, *qctxk;

        long *a, i, n;
        qadic_dense_t Delta_x, Delta_y, Delta, V, y;

        long M, *b, l, m;
        qadic_dense_ctx_struct *qctx;

        /* Number of steps */
        n = FLINT_CLOG2(N - 2 * 3) + 1;

        /* Precision needed at each step */
        a = flint_malloc(n * sizeof(long));

        i = 0;
        for (a[i = 0] = N; a[i] > 2 * 3 + 1; i++)
        {
            a[i + 1] = (a[i] + 1) / 2 + 3;
        }
        /* i = n - 1 */

        qadic_dense_init(Delta_x);
        qadic_dense_init(Delta_y);
        qadic_dense_init(Delta);
        qadic_dense_init(V);
        qadic_dense_init(y);

        qctxi = flint_malloc(sizeof(qadic_dense_ctx_struct));
        qctxk = flint_malloc(sizeof(qadic_dense_ctx_struct));

        n = (FLINT_CLOG2(N) + 1);
        b = flint_malloc((2 * n) * sizeof(long));
        qctx = flint_malloc((2 * n) * sizeof(qadic_dense_ctx_struct));

        /* Lifting */
        {
            qadic_dense_ctx_init_reduce(qctxi, ctx, a[i]);
            qadic_dense_ctx_init_reduce(qctxk, ctx, a[i]);

            qadic_dense_set(rop, op);
        }
        for (i--; i >= 0; i--)
        {
            /* Compute precisions for Artin--Schreier subroutine */
            M = (a[i] + 1)/2;

            l = 0;
            {
                b[l] = M;
            }
            for (; !(b[l] & 1L); l++)
            {
                b[l + 1] = b[l] >> 1;
            }
            for (; b[l] > 2; l += 2)
            {
                b[l + 1] = (b[l] >> 1) + 1;
                b[l + 2] = b[l] >> 1;
            }
            if (b[l] == 2L)
            {
                b[l + 1] = 1L;
                l++;
            }

            /* Contexts for Artin--Schreier subroutine */
            m = l + 1;
            l = 0;
            {
                qadic_dense_ctx_init_reduce(qctx + l, ctx, b[l]);
            }
            for (l++; l < m; l++)
            {
                qadic_dense_ctx_init_reduce(qctx + l, qctx + l - 1, b[l]);
            }

            /* Context for previous step */
            qctxj = qctxk;
            qctxk = qctxi;
            qctxi = qctxj;

            /* Context for current step */
            qadic_dense_ctx_clear(qctxi);
            qadic_dense_ctx_init_reduce(qctxi, ctx, a[i]);

            /* y */
            qadic_dense_frobenius_teichmuller(y, rop, qctxi);
            /* 1 */
            qadic_dense_one(Delta_x, qctxi);
            /* 1 + x */
            qadic_dense_add(Delta_x, Delta_x, rop, qctxi);
            /* x y */
            qadic_dense_mul(Delta_y, rop, y, qctxi);
            /* 4 x y */
            Delta_y->val += 2;
            padic_poly_reduce(Delta_y, &qctxi->pctx);

            /* V */
            /* 4 x y² */
            qadic_dense_mul(V, Delta_y, y, qctxi);
            /* (1 + x)² */
            qadic_dense_sqr(Delta, Delta_x, qctxi);
            /* 4 x y² - (1 + x)² */
            qadic_dense_sub(V, V, Delta, qctxi);
            /* 4 x y² - (1 + x)² / 2^N'*/
            V->val -= a[i + 1];

            /* Delta_x */
            /* 2 (1 + x) */
            Delta_x->val += 1;
            padic_poly_reduce(Delta_x, &qctxk->pctx);
            /* y */
            padic_poly_reduce(y, &qctxk->pctx);
            /* y² */
            qadic_dense_sqr(Delta, y, qctxk);
            /* 4 y² */
            Delta->val += 2;
            padic_poly_reduce(Delta, &qctxk->pctx);
            /* 4 y² - 2 (1 + x)² */
            qadic_dense_sub(Delta_x, Delta, Delta_x, qctxk);
            /* 4 y² - 2 (1 + x)² / 2^k */
            Delta_x->val -= 3;

            /* Delta_y */
            /* 8 x y */
            Delta_y->val += 1;
            padic_poly_reduce(Delta_y, &qctxk->pctx);
            /* 8 x y / 2^k */
            Delta_y->val -= 3;

            /* Delta */
            _qadic_dense_artin_schreier_root_ii(Delta, Delta_y, Delta_x, V, roots, m, b, qctx);
            Delta->val += M;

            /* rop */
            qadic_dense_add(rop, rop, Delta, qctxi);
            padic_poly_reduce(rop, &qctxi->pctx);

#if DEBUG >= 1
            /* y */
            qadic_dense_frobenius_teichmuller(y, rop, qctxi);
            /* 1 */
            qadic_dense_one(Delta_x, qctxi);
            /* 1 + x */
            qadic_dense_add(Delta_x, Delta_x, rop, qctxi);
            /* x y */
            qadic_dense_mul(Delta_y, rop, y, qctxi);
            /* 4 x y */
            Delta_y->val += 2;
            padic_poly_reduce(Delta_y, &qctxi->pctx);

            /* V */
            /* 4 x y² */
            qadic_dense_mul(V, Delta_y, y, qctxi);
            /* (1 + x)² */
            qadic_dense_sqr(Delta, Delta_x, qctxi);
            /* 4 x y² - (1 + x)² */
            qadic_dense_sub(V, V, Delta, qctxi);
            padic_poly_reduce(V, &qctxi->pctx);
            printf("Lambda'[%ld] = ", a[i]), qadic_dense_print_pretty(V, qctxi) ,printf("\n");
#endif

            for (l = 0; l < m; l++)
                qadic_dense_ctx_clear(qctx + l);
        }

        qadic_dense_clear(Delta_x);
        qadic_dense_clear(Delta_y);
        qadic_dense_clear(Delta);
        qadic_dense_clear(V);
        qadic_dense_clear(y);

        qadic_dense_ctx_clear(qctxi);
        qadic_dense_ctx_clear(qctxk);
        flint_free(qctxi);
        flint_free(qctxk);

        flint_free(b);
        flint_free(qctx);
    }
    else
    {
        printf("ERROR (_padic_poly_teichmuller).  odd characteristic not implemented.\n");
        abort();
    }
}

/* Newton lift for the modified modular polynomial \Lambda'_2 */
/* \Lambda'_2 = 4 X Y^2 - (1 + X)^2 */
/* k = 3 */
void qadic_dense_gen_newton_lift_ii_invlambda_2(qadic_dense_t rop, const qadic_dense_t op,
                                      const qadic_dense_struct *roots,
                                      const qadic_dense_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (N < 2 * 3 + 1 + 1)
    {
        qadic_dense_set(rop, op);
    }
    else if (*(&ctx->pctx)->p == 2L)
    {
        _qadic_dense_gen_newton_lift_ii_invlambda_2(rop, op, roots, N, ctx);
    }
    else
    {
        printf("ERROR (_padic_poly_teichmuller).  odd characteristic not implemented.\n");
        abort();
    }
}

void _qadic_dense_traces_init(padic_struct *traces, qadic_dense_ctx_t ctx)
{
    const long d = qadic_dense_ctx_degree(ctx);

    long i, j;
    padic_t m, n;

    padic_init(m, &ctx->pctx);
    padic_init(n, &ctx->pctx);

    i = 0;
    {
        padic_set_ui(traces + 0, d, &ctx->pctx);
    }
    for (i++ ; i < d; i++)
    {
        padic_set_si(n, -i, &ctx->pctx);
        padic_poly_get_coeff_padic(traces + i, ctx->mod, d - i, &ctx->pctx);
        padic_mul(traces + i, traces + i, n, &ctx->pctx);
        for (j = 1; j < i; j++)
        {
            padic_poly_get_coeff_padic(n, ctx->mod, d - j, &ctx->pctx);
            padic_mul(n, traces + i - j, n, &ctx->pctx);
            padic_sub(traces + i, traces + i, n, &ctx->pctx);
        }
    }

    padic_clear(m, &ctx->pctx);
    padic_clear(n, &ctx->pctx);
}

void qadic_dense_traces_init(padic_struct **traces, qadic_dense_ctx_t ctx)
{
    const long d = qadic_dense_ctx_degree(ctx);

    long i;

    *traces = flint_malloc(d * sizeof(padic_struct));

    for (i = 0; i < d; i++)
    {
        padic_init(*traces + i, &ctx->pctx);
    }

    _qadic_dense_traces_init(*traces, ctx);
}

void qadic_dense_trace_char_2(padic_t rop, const qadic_dense_t op, const padic_struct *traces, const qadic_dense_ctx_t ctx)
{
    const long d = qadic_dense_ctx_degree(ctx);

    long i;

    padic_t a;

    padic_init(a, &ctx->pctx);

    padic_zero(rop);

    for (i = 0; i < d; i++)
    {
        padic_poly_get_coeff_padic(a, op, i, &ctx->pctx);
        padic_mul(a, a, traces + i, &ctx->pctx);
        padic_add(rop, rop, a, &ctx->pctx);
    }

    padic_clear(a, &ctx->pctx);
}

void qadic_dense_norm_char_2(padic_t rop, const qadic_dense_t op, const padic_struct *traces, const qadic_dense_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;
    const long s = floor(sqrt(N)) / 2;

    long i, v;

    clock_t time;

    const fmpz_t two = {2L};
    fmpz_t pow;
    padic_t ptwo, r, rinv;

    padic_ctx_t pctx2, pctx3;
    padic_poly_t modulus;
    qadic_dense_ctx_t ctx2, ctx3;

    qadic_dense_t gamma, gamma2, t, w, z;

    fmpz_init(pow);
    fmpz_pow_ui(pow, two, s);

    padic_init(ptwo, &ctx->pctx);
    padic_set_ui(ptwo, 2, &ctx->pctx);
    padic_init(r, &ctx->pctx);
    padic_init(rinv, &ctx->pctx);

    padic_ctx_init(pctx2, two, N + s, (&ctx->pctx)->mode);
    padic_ctx_init(pctx3, two, N + s - 1, (&ctx->pctx)->mode);
    padic_poly_init(modulus);
    qadic_dense_ctx_modulus(modulus, ctx);
    qadic_dense_ctx_init(ctx2, modulus, pctx2, ctx->var);
    qadic_dense_ctx_init(ctx3, modulus, pctx3, ctx->var);

    qadic_dense_init(gamma);
    qadic_dense_init(gamma2);
    qadic_dense_init(t);
    qadic_dense_init(w);
    qadic_dense_init(z);

    /* z = (a^(2^s) - 1) / 2 */
    time = clock();
    qadic_dense_pow(z, op, pow, ctx2);
    qadic_dense_one(w, ctx2);
    qadic_dense_sub(z, z, w, ctx2);
    z->val -= 1;
    v = z->val;
    printf("pow took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);

    /* gamma = z / (1 + z)*/
    qadic_dense_add(gamma, z, w, ctx3);
    time = clock();
    qadic_dense_inv(gamma, gamma, ctx3);
    printf("inv took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
    qadic_dense_mul(gamma, gamma, z, ctx3);
    /* gamma2 = gamma^2*/
    qadic_dense_pow(gamma2, gamma, two, ctx3);

    /* w = 2^(-s) log(a) */
    time = clock();
    padic_one(r, pctx3);
    padic_poly_fit_length(w, 2*(gamma->length)*(N + s));
    padic_poly_zero(w);
    for (i = v; i < N + s - 1 - 2*v; i += 2*v)
    {
        padic_inv(rinv, r, pctx3);
        padic_poly_scalar_mul_padic(t, gamma, rinv, pctx3);
        padic_poly_reduce(t, pctx3);
        padic_poly_add(w, w, t, pctx3);
        qadic_dense_mul(gamma, gamma, gamma2, ctx3);
        padic_add(r, r, ptwo, pctx3);
    }
    {
        padic_inv(rinv, r, pctx3);
        padic_poly_scalar_mul_padic(t, gamma, rinv, pctx3);
        padic_poly_reduce(t, pctx3);
        padic_poly_add(w, w, t, pctx3);
    }
    w->val -= s - 1;
    printf("log took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);

    /* tr(w) */
    /* use precomputations? */
    /* qadic_dense_trace(rop, w, ctx); */
    time = clock();
    qadic_dense_trace_char_2(rop, w, traces, ctx);
    printf("trace took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);

    /* exp(tr(w)) */
    time = clock();
    padic_exp(rop, rop, &ctx->pctx);
    printf("exp took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
    qadic_dense_clear(gamma);
    qadic_dense_clear(gamma2);
    qadic_dense_clear(t);
    qadic_dense_clear(w);
    qadic_dense_clear(z);

    padic_ctx_clear(pctx2);
    padic_ctx_clear(pctx3);
    padic_poly_clear(modulus);
    qadic_dense_ctx_clear(ctx2);
    qadic_dense_ctx_clear(ctx3);

    padic_clear(ptwo, &ctx->pctx);
    padic_clear(r, &ctx->pctx);
    padic_clear(rinv, &ctx->pctx);

    fmpz_clear(pow);
}

int harley_gamma_2(long d, long N)
{
    clock_t time, total_time;
    long N_trace;
    fmpz_t p, pow, s, t2;
    padic_t t;
    padic_struct *traces;
    fmpz_poly_t f, g;
    padic_poly_t h;
    padic_ctx_t pctx, pctx_trace;
    qadic_dense_ctx_t qctx, qctx_one;
    qadic_dense_struct *roots;
    qadic_dense_t a, x, xinv;

    total_time = clock();

    fmpz_init(p);
    fmpz_set_ui(p, 2);

    N_trace = ((d + 1) >> 1) + 2;

    if (N == 0)
        N = N_trace;

    printf("d = %ld\n", d);
    printf("N = %ld\n", N);

    if (fmpz_poly_init_conway(f, p, d) && fmpz_poly_init_conway_mgm(f, p, d))
    {
        printf("Error.  Conway polynomial not found.\n");
        abort();
    }
    printf("f = ");
    fmpz_poly_print_pretty(f, "x");
    printf("\n");

    fmpz_poly_init2(g, d + 1);
    time = clock();
    _padic_poly_teichmuller(g, f, p, N);
    printf("Modulus took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
#if DEBUG >= 1
    printf("g = ");
    fmpz_poly_print_pretty(g, "x");
    printf("\n");
#endif

    time = clock();
    padic_ctx_init(pctx_trace, p, N_trace, PADIC_TERSE);
    padic_ctx_init(pctx, p, N, PADIC_TERSE);

    padic_poly_init(h);
    padic_poly_set_fmpz_poly(h, g, pctx);
    printf("p-adic context took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
#if DEBUG >= 1
    printf("h = ");
    padic_poly_print_pretty(h, "x", pctx);
    printf("\n");
#endif

    time = clock();
    qadic_dense_ctx_init(qctx, h, pctx, "t");
    qadic_dense_ctx_init_reduce(qctx_one, qctx, 1);
    printf("q-adic context took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
#if DEBUG >= 1
    printf("invmod = ");
    padic_poly_print_pretty(qctx->invmod, "x", pctx);
    printf("\n");
#endif

    time = clock();
    qadic_dense_proot_basis_init(&roots, qctx_one);
    printf("p-th root took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
#if DEBUG >= 1
    printf("roots[1] = ");
    qadic_dense_print_pretty(roots + 1, qctx);
    printf("\n");
#endif
    /*
    time = clock();
    qadic_dense_frobenius_teichmuller_precomp_init_char_2(&frobs, qctx);
    printf("Frobenius took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
#if DEBUG >= 1
    printf("frobs[1] = ");
    qadic_dense_print_pretty(frobs + 1, qctx);
    printf("\n");
#endif
    time = clock();
    qadic_dense_frobenius_teichmuller_precomp_reduce_char_2(&frobs, qctx);
    printf("Frobenius reduction took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
    */

    padic_poly_init2(a, d - 1);
    _padic_poly_set_length(a, d - 1);
    fmpz_set_ui(a->coeffs + d - 2, 1);
#if DEBUG >= 1
    printf("a = ");
    qadic_dense_print_pretty(a, qctx);
    printf("\n");
#endif
    time = clock();
    qadic_dense_init(x);
    qadic_dense_gen_newton_lift_ii_gamma_2(x, a, roots, qctx);
#if DEBUG >= 1
    printf("x = ");
    qadic_dense_print_pretty(x, qctx);
    printf("\n");
#endif
    printf("Lift took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
    abort();

    time = clock();
    qadic_dense_init(xinv);
    qadic_dense_set(xinv, x);
    /* groumpf, that means, we can gain 2 bits above... */
    xinv->val += 2;
    padic_poly_reduce(xinv, &qctx->pctx);
    qadic_dense_one(x, qctx);
    qadic_dense_add(xinv, xinv, x, qctx);
    qadic_dense_inv(xinv, xinv, qctx);
#if DEBUG >= 1
    printf("xinv = ");
    qadic_dense_print_pretty(xinv, qctx);
    printf("\n");
#endif
    printf("Inverse took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);

    /* Go back to sparse modulus from here... */
    time = clock();
    qadic_dense_traces_init(&traces, qctx);
    printf("Traces took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
    time = clock();
    padic_init(t, pctx);
    qadic_dense_norm_char_2(t, xinv, traces, qctx);
    /* qadic_dense_norm(t, xinv, qctx); */
    printf("norm = ");
    padic_print(t, pctx);
    printf("\n");
    printf("Norm took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);

    padic_reduce(t, pctx_trace);
    fmpz_init(t2);
    padic_get_fmpz(t2, t, pctx_trace);
    fmpz_init(s);
    fmpz_pow_ui(s, t2, 2);
    fmpz_init(pow);
    fmpz_pow_ui(pow, p, d + 1);

    if (fmpz_cmp(s, pow) > 0)
    {
        fmpz_pow_ui(pow, p, N_trace);
        fmpz_sub(t2, t2, pow);
    }
    printf("trace = ");
    fmpz_print(t2);
    printf("\n");

    fmpz_pow_ui(pow, p, d);
    fmpz_add_ui(pow, pow, 1);
    fmpz_sub(t2, pow, t2);
    printf("card = ");
    fmpz_print(t2);
    printf("\n");
    printf("Total time: %f seconds.\n", ((double) (clock() - total_time)) / CLOCKS_PER_SEC);

    padic_clear(t, pctx);

    qadic_dense_clear(a);
    qadic_dense_clear(x);
    qadic_dense_clear(xinv);

    qadic_dense_ctx_clear(qctx);
    padic_ctx_clear(pctx);
    padic_ctx_clear(pctx_trace);

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);

    fmpz_clear(p);
    fmpz_clear(pow);
    fmpz_clear(s);
    fmpz_clear(t2);

    return 0;
}

void qadic_dense_invsqrt_lambda(qadic_dense_t rop, const qadic_dense_t lambda, const qadic_dense_t a, const qadic_dense_ctx_t ctx)
{
    long i;
    qadic_dense_t one;
    qadic_dense_t x;

    qadic_dense_init(one);
    qadic_dense_one(one, ctx);
    qadic_dense_init(x);

    /* N = 3 */
    qadic_dense_set(rop, a);
    rop->val += 2;
    qadic_dense_sub(rop, one, rop, ctx);

    /* N = 5, 7 */
    for (i = 0; i < 2; i++)
    {
        qadic_dense_mul(x, rop, rop, ctx);
        qadic_dense_mul(x, x, lambda, ctx);
        qadic_dense_sub(x, one, x, ctx);
        x->val -= 1;
        qadic_dense_mul(x, x, rop, ctx);
        qadic_dense_add(rop, rop, x, ctx);
    }

    qadic_dense_clear(one);
    qadic_dense_clear(x);
}

int harley_lambda_2(long d, long N)
{
    clock_t time, total_time;
    long i;
    long N_trace;
    fmpz_t p, pow, s, t2;
    padic_t t;
    padic_struct *traces;
    fmpz_poly_t f, g;
    padic_poly_t h;
    padic_ctx_t pctx, pctx_trace;
    qadic_dense_ctx_t qctx, qctxk, qctx_one;
    qadic_dense_struct *roots;
    qadic_dense_t one, a, x, xinv;

    total_time = clock();

    fmpz_init(p);
    fmpz_set_ui(p, 2);

    N_trace = ((d + 1) >> 1) + 2;

    if (N == 0)
        N = N_trace + 4;

    printf("d = %ld\n", d);
    printf("N = %ld\n", N);

    if (fmpz_poly_init_conway(f, p, d) && fmpz_poly_init_conway_mgm(f, p, d))
    {
        printf("Error.  Conway polynomial not found.\n");
        abort();
    }
    printf("f = ");
    fmpz_poly_print_pretty(f, "x");
    printf("\n");

    fmpz_poly_init2(g, d + 1);
    time = clock();
    _padic_poly_teichmuller(g, f, p, N);
    printf("Modulus took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
#if DEBUG >= 1
    printf("g = ");
    fmpz_poly_print_pretty(g, "x");
    printf("\n");
#endif

    time = clock();
    padic_ctx_init(pctx_trace, p, N_trace, PADIC_TERSE);
    padic_ctx_init(pctx, p, N, PADIC_TERSE);

    padic_poly_init(h);
    padic_poly_set_fmpz_poly(h, g, pctx);
    printf("p-adic context took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
#if DEBUG >= 1
    printf("h = ");
    padic_poly_print_pretty(h, "x", pctx);
    printf("\n");
#endif

    time = clock();
    qadic_dense_ctx_init(qctx, h, pctx, "t");
    qadic_dense_ctx_init_reduce(qctx_one, qctx, 1);
    qadic_dense_ctx_init_reduce(qctxk, qctx, 2 * 3 + 1);
    qadic_dense_init(one);
    qadic_dense_one(one, qctxk);
    printf("q-adic context took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
#if DEBUG >= 1
    printf("invmod = ");
    padic_poly_print_pretty(qctx->invmod, "x", pctx);
    printf("\n");
#endif

    time = clock();
    qadic_dense_proot_basis_init(&roots, qctx_one);
    printf("p-th root took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
#if DEBUG >= 1
    printf("roots[1] = ");
    qadic_dense_print_pretty(roots + 1, qctx);
    printf("\n");
#endif
    /*
    time = clock();
    qadic_dense_frobenius_teichmuller_precomp_init_char_2(&frobs, qctx);
    printf("Frobenius took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
#if DEBUG >= 1
    printf("frobs[1] = ");
    qadic_dense_print_pretty(frobs + 1, qctx);
    printf("\n");
#endif
    time = clock();
    qadic_dense_frobenius_teichmuller_precomp_reduce_char_2(&frobs, qctx);
    printf("Frobenius reduction took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
    */

    padic_poly_init2(a, d - 1);
    _padic_poly_set_length(a, d - 1);
    fmpz_set_ui(a->coeffs + d - 2, 1);
#if DEBUG >= 1
    printf("a = ");
    qadic_dense_print_pretty(a, qctx);
    printf("\n");
#endif

    time = clock();
    /* Compute initial value to sufficient precision */
    qadic_dense_init(x);
    qadic_dense_init(xinv);
    /* lambda_1 */
    qadic_dense_set(x, a);
    x->val += 3;
    qadic_dense_add(x, x, one, qctxk);
    /* lambda_2,3,4 */
    for (i = 2; i < 5; i++)
    {
        qadic_dense_sub(a, x, one, qctxk);
        a->val -= 3;
        qadic_dense_invsqrt_lambda(xinv, x, a, qctxk);
        qadic_dense_add(x, x, one, qctxk);
        x->val -= 1;
        qadic_dense_mul(x, x, xinv, qctxk);
        qadic_dense_inv(x, x, qctxk);
    }
    /* Lift it */
    qadic_dense_gen_newton_lift_ii_lambda_2(xinv, x, roots, qctx);
#if DEBUG >= 1
    printf("x = ");
    qadic_dense_print_pretty(xinv, qctx);
    printf("\n");
#endif
    printf("Lift took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);

    time = clock();
    qadic_dense_add(xinv, xinv, one, qctx);
    xinv->val -= 1;
    qadic_dense_inv(xinv, xinv, qctx);
#if DEBUG >= 1
    printf("xinv = ");
    qadic_dense_print_pretty(xinv, qctx);
    printf("\n");
#endif
    printf("Inverse took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);

    /* Go back to sparse modulus from here... */
    time = clock();
    qadic_dense_traces_init(&traces, qctx);
    printf("Traces took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
    time = clock();
    padic_init(t, pctx);
    qadic_dense_norm_char_2(t, xinv, traces, qctx);
    /* qadic_dense_norm(t, xinv, qctx); */
    printf("norm = ");
    padic_print(t, pctx);
    printf("\n");
    printf("Norm took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);

    padic_reduce(t, pctx_trace);
    fmpz_init(t2);
    padic_get_fmpz(t2, t, pctx_trace);
    fmpz_init(s);
    fmpz_pow_ui(s, t2, 2);
    fmpz_init(pow);
    fmpz_pow_ui(pow, p, d + 1);

    if (fmpz_cmp(s, pow) > 0)
    {
        fmpz_pow_ui(pow, p, N_trace);
        fmpz_sub(t2, t2, pow);
    }
    printf("trace = ");
    fmpz_print(t2);
    printf("\n");

    fmpz_pow_ui(pow, p, d);
    fmpz_add_ui(pow, pow, 1);
    fmpz_sub(t2, pow, t2);
    printf("card = ");
    fmpz_print(t2);
    printf("\n");
    printf("Total time: %f seconds.\n", ((double) (clock() - total_time)) / CLOCKS_PER_SEC);

    padic_clear(t, pctx);

    qadic_dense_clear(one);
    qadic_dense_clear(a);
    qadic_dense_clear(x);
    qadic_dense_clear(xinv);

    qadic_dense_ctx_clear(qctx);
    padic_ctx_clear(pctx);
    padic_ctx_clear(pctx_trace);

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);

    fmpz_clear(p);
    fmpz_clear(pow);
    fmpz_clear(s);
    fmpz_clear(t2);

    return 0;
}

int harley_invlambda_2(long d, long N)
{
    clock_t time, total_time;
    long i;
    long N_trace;
    fmpz_t p, pow, s, t2;
    padic_t t;
    padic_struct *traces;
    fmpz_poly_t f, g;
    padic_poly_t h;
    padic_ctx_t pctx, pctx_trace;
    qadic_dense_ctx_t qctx, qctxk, qctx_one;
    qadic_dense_struct *roots;
    qadic_dense_t one, a, x, xinv;

    total_time = clock();

    fmpz_init(p);
    fmpz_set_ui(p, 2);

    N_trace = ((d + 1) >> 1) + 2;

    if (N == 0)
        N = N_trace + 4;

    printf("d = %ld\n", d);
    printf("N = %ld\n", N);

    if (fmpz_poly_init_conway(f, p, d) && fmpz_poly_init_conway_mgm(f, p, d))
    {
        printf("Error.  Conway polynomial not found.\n");
        abort();
    }
    printf("f = ");
    fmpz_poly_print_pretty(f, "x");
    printf("\n");

    fmpz_poly_init2(g, d + 1);
    time = clock();
    _padic_poly_teichmuller(g, f, p, N);
    printf("Modulus took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
#if DEBUG >= 1
    printf("g = ");
    fmpz_poly_print_pretty(g, "x");
    printf("\n");
#endif

    padic_ctx_init(pctx_trace, p, N_trace, PADIC_TERSE);
    padic_ctx_init(pctx, p, N, PADIC_TERSE);

    padic_poly_init(h);
    padic_poly_set_fmpz_poly(h, g, pctx);
#if DEBUG >= 1
    printf("h = ");
    padic_poly_print_pretty(h, "x", pctx);
    printf("\n");
#endif

    qadic_dense_ctx_init(qctx, h, pctx, "t");
    qadic_dense_ctx_init_reduce(qctx_one, qctx, 1);
    qadic_dense_ctx_init_reduce(qctxk, qctx, 2 * 3 + 1);
    qadic_dense_init(one);
    qadic_dense_one(one, qctxk);
#if DEBUG >= 1
    printf("invmod = ");
    padic_poly_print_pretty(qctx->invmod, "x", pctx);
    printf("\n");
#endif

    time = clock();
    qadic_dense_proot_basis_init(&roots, qctx_one);
    printf("p-th root took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
#if DEBUG >= 1
    printf("roots[1] = ");
    qadic_dense_print_pretty(roots + 1, qctx);
    printf("\n");
#endif
    /*
    time = clock();
    qadic_dense_frobenius_teichmuller_precomp_init_char_2(&frobs, qctx);
    printf("Frobenius took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
#if DEBUG >= 1
    printf("frobs[1] = ");
    qadic_dense_print_pretty(frobs + 1, qctx);
    printf("\n");
#endif
    time = clock();
    qadic_dense_frobenius_teichmuller_precomp_reduce_char_2(&frobs, qctx);
    printf("Frobenius reduction took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
    */

    padic_poly_init2(a, d - 1);
    _padic_poly_set_length(a, d - 1);
    fmpz_set_ui(a->coeffs + d - 2, 1);
#if DEBUG >= 1
    printf("a = ");
    qadic_dense_print_pretty(a, qctx);
    printf("\n");
#endif

    time = clock();
    /* Compute initial value to sufficient precision */
    qadic_dense_init(x);
    qadic_dense_init(xinv);
    /* lambda_1 */
    qadic_dense_set(x, a);
    x->val += 3;
    qadic_dense_add(x, x, one, qctxk);
    /* lambda_2,3,4 */
    for (i = 2; i < 5; i++)
    {
        qadic_dense_sub(a, x, one, qctxk);
        a->val -= 3;
        qadic_dense_invsqrt_lambda(xinv, x, a, qctxk);
        qadic_dense_add(x, x, one, qctxk);
        x->val -= 1;
        qadic_dense_mul(x, x, xinv, qctxk);
    }
    /* Lift it */
    qadic_dense_gen_newton_lift_ii_invlambda_2(xinv, x, roots, qctx);
#if DEBUG >= 1
    printf("xinv = ");
    qadic_dense_print_pretty(xinv, qctx);
    printf("\n");
#endif
    printf("Lift took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);

    time = clock();
    qadic_dense_add(x, xinv, one, qctx);
    x->val -= 1;
    qadic_dense_inv(x, x, qctx);
    qadic_dense_mul(xinv, xinv, x, qctx);
#if DEBUG >= 1
    printf("xinv = ");
    qadic_dense_print_pretty(xinv, qctx);
    printf("\n");
#endif
    printf("Inverse took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);

    /* Go back to sparse modulus from here... */
    time = clock();
    qadic_dense_traces_init(&traces, qctx);
    printf("Traces took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
    time = clock();
    padic_init(t, pctx);
    qadic_dense_norm_char_2(t, xinv, traces, qctx);
    /* qadic_dense_norm(t, xinv, qctx); */
    printf("norm = ");
    padic_print(t, pctx);
    printf("\n");
    printf("Norm took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);

    padic_reduce(t, pctx_trace);
    fmpz_init(t2);
    padic_get_fmpz(t2, t, pctx_trace);
    fmpz_init(s);
    fmpz_pow_ui(s, t2, 2);
    fmpz_init(pow);
    fmpz_pow_ui(pow, p, d + 1);

    if (fmpz_cmp(s, pow) > 0)
    {
        fmpz_pow_ui(pow, p, N_trace);
        fmpz_sub(t2, t2, pow);
    }
    printf("trace = ");
    fmpz_print(t2);
    printf("\n");

    fmpz_pow_ui(pow, p, d);
    fmpz_add_ui(pow, pow, 1);
    fmpz_sub(t2, pow, t2);
    printf("card = ");
    fmpz_print(t2);
    printf("\n");
    printf("Total time: %f seconds.\n", ((double) (clock() - total_time)) / CLOCKS_PER_SEC);

    padic_clear(t, pctx);

    qadic_dense_clear(one);
    qadic_dense_clear(a);
    qadic_dense_clear(x);
    qadic_dense_clear(xinv);

    qadic_dense_ctx_clear(qctx);
    padic_ctx_clear(pctx);
    padic_ctx_clear(pctx_trace);

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);

    fmpz_clear(p);
    fmpz_clear(pow);
    fmpz_clear(s);
    fmpz_clear(t2);

    return 0;
}

int main(int argc, char *argv[])
{
    long d, N;

    if (argc < 2)
    {
        printf("Usage: points d [N]\n");
        return 0;
    }

    d = atol(argv[1]);

    if (argc > 2)
        N = atol(argv[2]);
    else
        N = 0;

    return harley_gamma_2(d, N);
    /*return (harley_gamma_2(d, N) || harley_lambda_2(d, N) || harley_invlambda_2(d, N));*/
}
