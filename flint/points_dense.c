#undef ulong
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define ulong unsigned long

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

#define MGM_LENGTH 31
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

    ctx->var = flint_malloc(strlen(var));
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

    padic_poly_init2(rop->invmod, d + 1);
    _fmpz_vec_scalar_mod_fmpz(rop->invmod->coeffs, op->invmod->coeffs, d + 1, pow);
    _padic_poly_set_length(rop->invmod, d + 1);
    _padic_poly_normalise(rop->invmod);

    rop->var = flint_malloc(strlen(op->var));
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

    padic_poly_init2(rop->invmod, d + 1);
    _fmpz_vec_scalar_fdiv_r_2exp(rop->invmod->coeffs, op->invmod->coeffs, d + 1, N);
    _padic_poly_set_length(rop->invmod, d + 1);
    _padic_poly_normalise(rop->invmod);

    rop->var = flint_malloc(strlen(op->var));
    strcpy(rop->var, op->var);

    return;
}

/* Powers of p should be shared */
void _padic_poly_teichmuller_inc_recursive(fmpz_poly_t delta, const fmpz_poly_t f0,
                                           const fmpz_poly_t f1, const fmpz_poly_t V,
                                           long d, long m, const long *b, const fmpz *pow)
{
    if (m == 1)
    {
        fmpz_poly_neg(delta, V);
        _fmpz_vec_scalar_mod_fmpz(delta->coeffs, delta->coeffs, delta->length, pow);
        _fmpz_poly_normalise(delta);
    }
    else if (*(pow + m  - 1) == 2L)
    {
        long *a, h, i, j, k, l, l0, l1, ld, n;
        const fmpz *powi;
        fmpz_poly_t delta0, df0, delta1, df1, W, Delta;

        /* Number of steps */
        n = FLINT_CLOG2(b[0]) + 1;

        /* Precision needed at each step */
        a = flint_malloc(n * sizeof(long));

        i = 0;
        for (a[i = 0] = b[0]; a[i] > 1; i++)
        {
            a[i + 1] = (a[i] + 1) / 2;
        }

        /* Degree and length */
        l = d + 1;
        l0 = (l + 1) >> 1;
        l1 = (l >> 1);

        fmpz_poly_init2(delta0, l0);
        fmpz_poly_init2(df0, l);
        fmpz_poly_init2(delta1, l1);
        fmpz_poly_init2(df1, l);
        fmpz_poly_init2(W, l);
        fmpz_poly_init2(Delta, l);

        /* Lifting */
        i = n - 1;
        j = m - 1;
        {
            powi = pow + j;

            fmpz_poly_set(delta, V);

            ld = delta->length;

            _fmpz_vec_scalar_mod_fmpz(delta->coeffs, delta->coeffs, ld, powi);
            _fmpz_poly_normalise(delta);
        }
        for (i--; i >= 0; i--)
        {
            h = j;
            j--;
            powi--;
            if (b[j] != a[i])
            {
                j--;
                powi--;
            }

            ld = delta->length;

            /* delta_0, delta_1, and products */
            _fmpz_poly_set_length(delta0, (ld + 1) >> 1);
            _fmpz_poly_set_length(delta1, (ld >> 1));

            for (k = 0; k < ld; k++)
                fmpz_set(((k % 2)?delta1:delta0)->coeffs + (k >> 1), delta->coeffs + k);

            _fmpz_poly_normalise(delta0);
            _fmpz_poly_normalise(delta1);

            _fmpz_vec_scalar_mod_fmpz(df0->coeffs, f0->coeffs, f0->length, powi);
            _fmpz_vec_scalar_mod_fmpz(df1->coeffs, f1->coeffs, f1->length, powi);
            _fmpz_poly_set_length(df0, f0->length);
            _fmpz_poly_set_length(df1, f1->length);
            _fmpz_poly_normalise(df0);
            _fmpz_poly_normalise(df1);

            _fmpz_mod_poly_mul(df0->coeffs, delta0->coeffs, delta0->length, f0->coeffs, f0->length, powi);
            _fmpz_mod_poly_mul(df1->coeffs, delta1->coeffs, delta1->length, f1->coeffs, f1->length, powi);

            _fmpz_poly_set_length(df0, delta0->length + f0->length - 1);
            _fmpz_poly_set_length(df1, delta1->length + f1->length);

            _fmpz_poly_shift_left(df1->coeffs, df1->coeffs, df1->length - 1, 1);

            _fmpz_mod_poly_sub(df0->coeffs, df0->coeffs, df0->length, df1->coeffs, df1->length, powi);

            _fmpz_poly_set_length(df0, FLINT_MAX(df0->length, df1->length));
            _fmpz_poly_normalise(df0);

            _fmpz_vec_scalar_mul_2exp(df0->coeffs, df0->coeffs, df0->length, 1);
            /* Computed everything above with one superfluous bit... */
            _fmpz_vec_scalar_mod_fmpz(df0->coeffs, df0->coeffs, df0->length, powi);
            _fmpz_poly_normalise(df0);

            /* W */
            _fmpz_mod_poly_add(W->coeffs, delta->coeffs, ld, V->coeffs, V->length, powi);
            _fmpz_poly_set_length(W, FLINT_MAX(ld, V->length));
            _fmpz_poly_normalise(W);

            _fmpz_mod_poly_sub(W->coeffs, W->coeffs, W->length, df0->coeffs, df0->length, powi);
            _fmpz_poly_set_length(W, FLINT_MAX(W->length, df0->length));
            _fmpz_poly_normalise(W);

            _fmpz_vec_scalar_fdiv_q_2exp(W->coeffs, W->coeffs, W->length, a[i + 1]);

            /* Delta */
            if (b[h] == a[i] - a[i + 1])
                _padic_poly_teichmuller_inc_recursive(Delta, f0, f1, W, d, m - h, b + h, pow + h);
            else
                _padic_poly_teichmuller_inc_recursive(Delta, f0, f1, W, d, m - h - 1, b + h + 1, pow + h + 1);
            _fmpz_vec_scalar_mul_2exp(Delta->coeffs, Delta->coeffs, Delta->length, a[i + 1]);
            /* Check precision here */
            _fmpz_vec_scalar_mod_fmpz(Delta->coeffs, Delta->coeffs, Delta->length, powi);
            _fmpz_poly_normalise(Delta);

            /* update delta */
            _fmpz_mod_poly_add(delta->coeffs, delta->coeffs, ld, Delta->coeffs, Delta->length, powi);
            _fmpz_poly_set_length(delta, FLINT_MAX(ld, Delta->length));
            _fmpz_poly_normalise(delta);
        }

        fmpz_poly_clear(delta0);
        fmpz_poly_clear(df0);
        fmpz_poly_clear(delta1);
        fmpz_poly_clear(df1);
        fmpz_poly_clear(W);
        fmpz_poly_clear(Delta);
    }
}

void _padic_poly_teichmuller_inc(fmpz_poly_t delta, const fmpz_poly_t f0,
                                 const fmpz_poly_t f1, const fmpz_poly_t V,
                                 long d, const fmpz_t p, long N)
{
    if (N == 1)
    {
        fmpz_poly_neg(delta, V);
        _fmpz_vec_scalar_mod_fmpz(delta->coeffs, delta->coeffs, delta->length, p);
        _fmpz_poly_normalise(delta);
    }
    else if (*(p) == 2L)
    {
        long *b, i, m, n;
        fmpz *pow;

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

        m = i + 1;
        pow = _fmpz_vec_init(m);

        /* Compute powers of 2 needed later */
        {
            fmpz_set_ui(pow + i, 2L);
        }
        for (i--; i >= 0; i--)
        {
            fmpz_mul_2exp(pow + i, pow + m - 1, b[i] - 1);
        }

        /* Let's go! */
        _padic_poly_teichmuller_inc_recursive(delta, f0, f1, V, d, m, b, pow);

        /* Deallocation */
        _fmpz_vec_clear(pow, m);
        flint_free(b);
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
        /* One should be sufficient */
        fmpz *pow, *powi;
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
            a[i + 1] = (a[i] + 1) / 2;
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
        pow = _fmpz_vec_init(m);

        /* Compute powers of 2 needed later */
        {
            fmpz_set_ui(pow + i, 2L);
        }
        for (i--; i >= 0; i--)
        {
            fmpz_mul_2exp(pow + i, pow + m - 1, b[i] - 1);
        }

        /* Lifting */
        i = n - 1;
        j = m - 1;
        {
            powi = pow + j;

            fmpz_poly_set(g, f);
        }
        for (i--; i >= 0; i--)
        {
            h = j;
            j--;
            powi--;
            if (b[j] != a[i])
            {
                j--;
                powi--;
            }

            /* f_0, f_1, and squares */
            _fmpz_poly_set_length(f0, l0);
            _fmpz_poly_set_length(f1, l1);

            for (k = 0; k < l; k++)
                fmpz_set(((k % 2)?f1:f0)->coeffs + (k >> 1), g->coeffs + k);

            _fmpz_poly_normalise(f0);
            _fmpz_poly_normalise(f1);

            _fmpz_mod_poly_sqr(f0_sqr->coeffs, f0->coeffs, f0->length, powi);
            _fmpz_mod_poly_sqr(f1_sqr->coeffs, f1->coeffs, f1->length, powi);

            _fmpz_poly_set_length(f0_sqr, (f0->length << 1) - 1);
            _fmpz_poly_set_length(f1_sqr, (f1->length << 1));

            _fmpz_poly_shift_left(f1_sqr->coeffs, f1_sqr->coeffs, f1_sqr->length - 1, 1);

            /* V */
            _fmpz_poly_set_length(V, l);

            _fmpz_vec_set(V->coeffs, g->coeffs, l);

            _fmpz_mod_poly_sub(V->coeffs, V->coeffs, V->length, f0_sqr->coeffs, f0_sqr->length, powi);
            _fmpz_poly_normalise(V);

            _fmpz_mod_poly_add(V->coeffs, V->coeffs, V->length, f1_sqr->coeffs, f1_sqr->length, powi);
            _fmpz_poly_set_length(V, FLINT_MAX(V->length, f1_sqr->length));
            _fmpz_poly_normalise(V);

            _fmpz_vec_scalar_fdiv_q_2exp(V->coeffs, V->coeffs, V->length, a[i + 1]);

            /* delta */
            if (b[h] == a[i] - a[i + 1])
                _padic_poly_teichmuller_inc_recursive(delta, f0, f1, V, d, m - h, b + h, pow + h);
            else
                _padic_poly_teichmuller_inc_recursive(delta, f0, f1, V, d, m - h - 1, b + h + 1, pow + h + 1);

            /* Use vec method */
            _fmpz_vec_scalar_mul_2exp(delta->coeffs, delta->coeffs, delta->length, a[i + 1]);
            _fmpz_vec_scalar_mod_fmpz(delta->coeffs, delta->coeffs, delta->length, powi);
            _fmpz_poly_normalise(delta);

            /* update g */
            _fmpz_mod_poly_add(g->coeffs, g->coeffs, l, delta->coeffs, delta->length, powi);
        }

        if (!fmpz_is_one(fmpz_poly_lead(g)))
            _fmpz_mod_poly_neg(g->coeffs, g->coeffs, g->length, powi);

        fmpz_poly_clear(f0);
        fmpz_poly_clear(f0_sqr);
        fmpz_poly_clear(f1);
        fmpz_poly_clear(f1_sqr);
        fmpz_poly_clear(V);
        fmpz_poly_clear(delta);

        _fmpz_vec_clear(pow, m);

        flint_free(a);
        flint_free(b);
    }
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
    else
    {
        long i;
        *roots = flint_malloc(p * sizeof(qadic_dense_struct));
        for (i = 0; i < p; i++)
            qadic_dense_init(*roots + i);
        _qadic_dense_proot_basis_init(*roots, ctx);
    }
}

void _qadic_dense_proot(qadic_dense_t rop, const qadic_dense_t op, const qadic_dense_struct *roots,
                  const qadic_dense_ctx_t ctx)
{
    const long d = qadic_dense_ctx_degree(ctx);
    const long p = *(&ctx->pctx)->p;
    long j, k;
    qadic_dense_t t;

    qadic_dense_init(t);

    qadic_dense_zero(rop);

    for (j = 0; j < p; j++)
    {
        for (k = 0; p*k + j < FLINT_MIN(d, op->length); k++)
            ;
        padic_poly_fit_length(t, k);
        for (k = 0; p*k + j < FLINT_MIN(d, op->length); k++)
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

void qadic_dense_frobenius_teichmuller_precomp_init(qadic_dense_struct **frobs,
                                      const qadic_dense_ctx_t ctx)
{
    const long d = qadic_dense_ctx_degree(ctx);

    long i;
    qadic_dense_t t;

    padic_poly_init2(t, d);
    fmpz_set_ui(t->coeffs, 0);
    fmpz_set_ui(t->coeffs + 1, 1);
    _padic_poly_set_length(t, 2);
    qadic_dense_pow(t, t, (&ctx->pctx)->p, ctx);

    *frobs = flint_malloc(d * sizeof(qadic_dense_struct));

    i = 0;
    {
        qadic_dense_init(*frobs);
        qadic_dense_one(*frobs, ctx);
    }
    for (i++; i < d; i++)
    {
        qadic_dense_init(*frobs + i);
        qadic_dense_mul(*frobs + i, *frobs + i - 1, t, ctx);
#if DEBUG >= 2
        printf("frobs[%ld] = ", i);
        qadic_dense_print_pretty(*frobs + i, ctx);
        printf("\n");
#endif
    }
}

void qadic_dense_frobenius_teichmuller_precomp_init_char_2(qadic_dense_struct **frobs,
                                      const qadic_dense_ctx_t ctx)
{
    /* The leading coefficient of the modulus is one */
    const long d = qadic_dense_ctx_degree(ctx);

    long i;
    long e;
    padic_poly_t r, s, t;
    padic_t a;

    /* r = x^d % mod */
    padic_poly_init2(r, d);
    _fmpz_vec_set(r->coeffs, ctx->mod->coeffs, d);
    _padic_poly_set_length(r, d);
    _padic_poly_normalise(s);
    padic_poly_neg(r, r, &ctx->pctx);
#if DEBUG >= 3
    printf("r = ");
    qadic_dense_print_pretty(r, ctx);
    printf("\n");
#endif
    /* s = x^(d+1) % mod */
    padic_poly_init2(s, d);
    _fmpz_vec_scalar_mul_fmpz(s->coeffs, ctx->mod->coeffs, d, ctx->mod->coeffs + d - 1);
    _fmpz_vec_scalar_fdiv_r_2exp(s->coeffs, s->coeffs, d, (&ctx->pctx)->N);
    _fmpz_vec_sub(s->coeffs + 1, s->coeffs + 1, ctx->mod->coeffs, d - 1);
    _fmpz_vec_scalar_fdiv_r_2exp(s->coeffs, s->coeffs, d, (&ctx->pctx)->N);
    _padic_poly_set_length(s, d);
    _padic_poly_normalise(s);
#if DEBUG >= 3
    printf("s = ");
    qadic_dense_print_pretty(s, ctx);
    printf("\n");
#endif

    padic_poly_init2(t, d);
    padic_init(a, &ctx->pctx);

    *frobs = flint_malloc(d * sizeof(qadic_dense_struct));

    i = 0;
    {
        qadic_dense_init(*frobs);
        qadic_dense_one(*frobs, ctx);
    }
    for (i++; i < d; i++)
    {
        e = FLINT_MIN(padic_poly_degree(*frobs + i - 1) + 1, d - 2);
        padic_poly_init2(*frobs + i, e + 2);
        fmpz_zero((*frobs + i)->coeffs);
        fmpz_zero((*frobs + i)->coeffs + 1);
        _fmpz_vec_set((*frobs + i)->coeffs + 2, (*frobs + i - 1)->coeffs, e);
        _padic_poly_set_length(*frobs + i, e + 2);
        _padic_poly_normalise(*frobs + i);

        padic_poly_get_coeff_padic(a, *frobs + i - 1, d - 2, &ctx->pctx);
        padic_poly_scalar_mul_padic(t, r, a, &ctx->pctx);
        padic_poly_add(*frobs + i, *frobs + i, t, &ctx->pctx);

        padic_poly_get_coeff_padic(a, *frobs + i - 1, d - 1, &ctx->pctx);
        padic_poly_scalar_mul_padic(t, s, a, &ctx->pctx);
        padic_poly_add(*frobs + i, *frobs + i, t, &ctx->pctx);

#if DEBUG >= 3
        printf("frobs[%ld] = ", i);
        qadic_dense_print_pretty(*frobs + i, ctx);
        printf("\n");
#endif
    }

    padic_poly_clear(t);
    padic_clear(a, &ctx->pctx);
}

void qadic_dense_frobenius_teichmuller_precomp_reduce(qadic_dense_struct **frobs,
                                              const qadic_dense_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (N == 1)
    {
        return;
    }
    else if (*(&ctx->pctx)->p == 2L)
    {
        const long d = qadic_dense_ctx_degree(ctx);

        long b, i, j, l, n;
        fmpz_t pow;

        fmpz_init(pow);

        /* Number of steps */
        n = FLINT_CLOG2(N) + 1;

        /* Compute precisions needed later */
        *frobs = flint_realloc(*frobs, d * (2 * n) * sizeof(qadic_dense_struct));

        i = 0;
        {
            b = N;
        }
        for (i++; !(b & 1L); i++)
        {
            b >>= 1;
            fmpz_one(pow);
            fmpz_mul_2exp(pow, pow, b);
            for (j = 0; j < d; j++)
            {
                l = padic_poly_length(*frobs + d * (i - 1) + j);
                padic_poly_init2(*frobs + d * i + j, l);
                _padic_poly_set_length(*frobs + d * i + j, l);
                _fmpz_vec_scalar_mod_fmpz((*frobs + d * i + j)->coeffs, (*frobs + d * (i - 1) + j)->coeffs, l, pow);
                _padic_poly_normalise(*frobs + d * i + j);
            }
        }
        for (; b > 2;i++)
        {
            b = (b >> 1) + 1;
            fmpz_one(pow);
            fmpz_mul_2exp(pow, pow, b);
            for (j = 0; j < d; j++)
            {
                l = padic_poly_length(*frobs + d * (i - 1) + j);
                padic_poly_init2(*frobs + d * i + j, l);
                _padic_poly_set_length(*frobs + d * i + j, l);
                _fmpz_vec_scalar_mod_fmpz((*frobs + d * i + j)->coeffs, (*frobs + d * (i - 1) + j)->coeffs, l, pow);
                _padic_poly_normalise(*frobs + d * i + j);
            }
            i++;
            b -= 1;
            fmpz_one(pow);
            fmpz_mul_2exp(pow, pow, b);
            for (j = 0; j < d; j++)
            {
                l = padic_poly_length(*frobs + d * (i - 1) + j);
                padic_poly_init2(*frobs + d * i + j, l);
                _padic_poly_set_length(*frobs + d * i + j, l);
                _fmpz_vec_scalar_mod_fmpz((*frobs + d * i + j)->coeffs, (*frobs + d * (i - 1) + j)->coeffs, l, pow);
                _padic_poly_normalise(*frobs + d * i + j);
            }
        }
        if (b == 2L)
        {
            b = 1L;
            fmpz_one(pow);
            fmpz_mul_2exp(pow, pow, b);
            for (j = 0; j < d; j++)
            {
                l = padic_poly_length(*frobs + d * (i - 1) + j);
                padic_poly_init2(*frobs + d * i + j, l);
                _padic_poly_set_length(*frobs + d * i + j, l);
                _fmpz_vec_scalar_mod_fmpz((*frobs + d * i + j)->coeffs, (*frobs + d * (i - 1) + j)->coeffs, l, pow);
                _padic_poly_normalise(*frobs + d * i + j);
            }
        }

        fmpz_clear(pow);
    }
}

void qadic_dense_frobenius_teichmuller__precomp_reduce_char_2(qadic_dense_struct **frobs,
                                                     const qadic_dense_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (N == 1)
    {
        return;
    }
    else if (*(&ctx->pctx)->p == 2L)
    {
        const long d = qadic_dense_ctx_degree(ctx);

        long b, i, j, l, n;

        /* Number of steps */
        n = FLINT_CLOG2(N) + 1;

        /* Compute precisions needed later */
        *frobs = flint_realloc(*frobs, d * (2 * n) * sizeof(qadic_dense_struct));

        i = 0;
        {
            b = N;
        }
        for (i++; !(b & 1L); i++)
        {
            b >>= 1;
            for (j = 0; j < d; j++)
            {
                l = padic_poly_length(*frobs + d * (i - 1) + j);
                padic_poly_init2(*frobs + d * i + j, l);
                _padic_poly_set_length(*frobs + d * i + j, l);
                _fmpz_vec_scalar_fdiv_r_2exp((*frobs + d * i + j)->coeffs, (*frobs + d * (i - 1) + j)->coeffs, l, b);
                _padic_poly_normalise(*frobs + d * i + j);
            }
        }
        for (; b > 2;i++)
        {
            b = (b >> 1) + 1;
            for (j = 0; j < d; j++)
            {
                l = padic_poly_length(*frobs + d * (i - 1) + j);
                padic_poly_init2(*frobs + d * i + j, l);
                _padic_poly_set_length(*frobs + d * i + j, l);
                _fmpz_vec_scalar_fdiv_r_2exp((*frobs + d * i + j)->coeffs, (*frobs + d * (i - 1) + j)->coeffs, l, b);
                _padic_poly_normalise(*frobs + d * i + j);
            }
            i++;
            b -= 1;
            for (j = 0; j < d; j++)
            {
                l = padic_poly_length(*frobs + d * (i - 1) + j);
                padic_poly_init2(*frobs + d * i + j, l);
                _padic_poly_set_length(*frobs + d * i + j, l);
                _fmpz_vec_scalar_fdiv_r_2exp((*frobs + d * i + j)->coeffs, (*frobs + d * (i - 1) + j)->coeffs, l, b);
                _padic_poly_normalise(*frobs + d * i + j);
            }
        }
        if (b == 2L)
        {
            b = 1L;
            for (j = 0; j < d; j++)
            {
                l = padic_poly_length(*frobs + d * (i - 1) + j);
                padic_poly_init2(*frobs + d * i + j, l);
                _padic_poly_set_length(*frobs + d * i + j, l);
                _fmpz_vec_scalar_fdiv_r_2exp((*frobs + d * i + j)->coeffs, (*frobs + d * (i - 1) + j)->coeffs, l, b);
                _padic_poly_normalise(*frobs + d * i + j);
            }
        }
    }
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

void _padic_poly_add_2exp(fmpz *rop, long *val,
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

void _qadic_dense_frobenius_teichmuller_precomp(qadic_dense_t rop, const qadic_dense_t op,
                                 const qadic_dense_struct *frobs,
                                 const qadic_dense_ctx_t ctx)
{
    const long d = qadic_dense_ctx_degree(ctx);

    long i;
    qadic_dense_t t;
    padic_t c;

    padic_init(c, &ctx->pctx);
    qadic_dense_init(t);

    padic_poly_fit_length(rop, d);
    qadic_dense_zero(rop);

    for (i = 0; i < d; i++)
    {
        padic_poly_get_coeff_padic(c, op, i, &ctx->pctx);
        padic_poly_scalar_mul_padic(t, frobs + i, c, &ctx->pctx);

        /* qadic_dense_add(rop, rop, t, ctx); */
        _padic_poly_add_no_mod(rop->coeffs, &(rop->val), rop->coeffs, rop->val, rop->length,
           t->coeffs, t->val, t->length, &ctx->pctx);
        /* _padic_poly_add_2exp(rop->coeffs, &(rop->val), rop->coeffs, rop->val, rop->length,
           t->coeffs, t->val, t->length, &ctx->pctx); */
        _padic_poly_set_length(rop, FLINT_MAX(rop->length, t->length));
    }
    padic_poly_reduce(rop, &ctx->pctx);
}

void _qadic_dense_frobenius_teichmuller_char_2(qadic_dense_t rop, const qadic_dense_t op,
                                 const qadic_dense_ctx_t ctx)
{
    const long len = padic_poly_length(op);
    long i;

    padic_poly_fit_length(rop, 2*len - 1);
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

void _qadic_dense_frobenius_teichmuller(qadic_dense_t rop, const qadic_dense_t op,
                                 const qadic_dense_ctx_t ctx)
{
    const long len = padic_poly_length(op);
    const long p = *(&ctx->pctx)->p;
    long i;

    padic_poly_fit_length(rop, p*len - 1);
    _fmpz_vec_zero(rop->coeffs, p*len - 1);
    for (i = 0; i < len; i++)
    {
        fmpz_set(rop->coeffs + p*i, op->coeffs + i);
    }
    _padic_poly_set_length(rop, p*len - 1);
    rop->val = op->val;

    qadic_dense_reduce(rop, ctx);
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
        qadic_dense_inv(x, alpha, qctx);
        qadic_dense_mul(x, x, gamma, qctx);
        qadic_dense_neg(x, x, qctx);
        qadic_dense_proot(x, x, roots, qctx);
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
            qctxi = qctx + j;

            qadic_dense_set(alpha2, alpha);
            qadic_dense_set(beta2, beta);
            qadic_dense_set(gamma2, gamma);

            padic_poly_reduce(alpha2, &qctxi->pctx);
            padic_poly_reduce(beta2, &qctxi->pctx);
            padic_poly_reduce(gamma2, &qctxi->pctx);

            qadic_dense_inv(x, alpha2, qctxi);
            qadic_dense_mul(x, x, gamma2, qctxi);
            qadic_dense_neg(x, x, qctxi);
            qadic_dense_proot(x, x, roots, qctxi);
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

            qadic_dense_frobenius_teichmuller(gamma3, x, qctxi);
            qadic_dense_mul(gamma3, gamma3, alpha2, qctxi);
            qadic_dense_mul(Delta2, beta2, x, qctxi);
            qadic_dense_add(gamma3, gamma3, Delta2, qctxi);
            qadic_dense_add(gamma3, gamma3, gamma2, qctxi);
            gamma3->val -= a[i + 1];

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
        qadic_dense_inv(x, alpha, ctx);
        qadic_dense_mul(x, x, gamma, ctx);
        qadic_dense_neg(x, x, ctx);
        qadic_dense_proot(x, x, roots, ctx);
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

/* Newton lift for modular polynomials */
/* k = 0, M = N' */
void _qadic_dense_gen_newton_lift_ii_modular(qadic_dense_t rop, const qadic_dense_t op,
                                      const qadic_dense_struct *roots,
                                      long m, const long *b, const qadic_dense_ctx_struct *qctx)
{
    if (m == 1)
    {
        qadic_dense_set(rop, op);
    }
    else if (*(&(qctx + m - 1)->pctx)->p == 2L)
    {
        const qadic_dense_ctx_struct *qctxi, *qctxk;
        long two;
        long *a, i, j, k, n;
        qadic_dense_t Delta_x, Delta_y, Delta, V, y, xy;

        two = 2L;

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
            k = j;
            qctxk = qctxi;
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
            qadic_dense_pow(V, Delta, &two, qctxi);
            /* (x + 2 y + 8 x y)² + y */
            qadic_dense_add(V, V, y, qctxi);
            /* 4 x y */
            xy->val += 1;
            padic_poly_reduce(xy, &qctxi->pctx);
            /* (x + 2 y + 8 x y)² + y + 4 x y */
            qadic_dense_add(V, V, xy, qctxi);
            /* (x + 2 y + 8 x y)² + y + 4 x y / 2^N' */
            V->val -= a[i + 1];

            /* Delta_x */
            /* x + 2 y + 8 x y */
            padic_poly_reduce(Delta, &qctxk->pctx);
            /* 1 */
            qadic_dense_one(xy, qctxk);
            /* 8 y */
            qadic_dense_set(Delta_x, y);
            Delta_x->val += 3;
            padic_poly_reduce(Delta_x, &qctxk->pctx);
            /* 1 + 8 y */
            qadic_dense_add(xy, xy, Delta_x, qctxk);
            /* (x + 2 y + 8 x y) (1 + 8 y) */
            qadic_dense_mul(Delta_x, Delta, xy, qctxk);
            /* 2 (x + 2 y + 8 x y) (1 + 8 y) */
            Delta_x->val += 1;
            padic_poly_reduce(Delta_x, &qctxk->pctx);
            /* 4 y */
            y->val += 2;
            padic_poly_reduce(y, &qctxk->pctx);
            /* 2 (x + 2 y + 8 x y) (1 + 8 y) + 4 y*/
            qadic_dense_add(Delta_x, Delta_x, y, qctxk);

            /* Delta_y */
            /* (x + 2 y + 8 x y) */
            Delta->val += 2;
            padic_poly_reduce(Delta, &qctxk->pctx);
            /* 1 + 2 (x + 2 y + 8 x y) */
            qadic_dense_one(xy, qctxk);
            qadic_dense_add(Delta_y, Delta, xy, qctxk);
            /* 4 x */
            qadic_dense_set(y, rop);
            y->val += 2;
            padic_poly_reduce(y, &qctxk->pctx);
            /* 1 + 4 x */
            qadic_dense_add(xy, xy, y, qctxk);
            /* (1 + 2 (x + 2 y + 8 x y)) (1 + 4 x) */
            qadic_dense_mul(Delta_y, Delta_y, xy, qctxk);

            /* Delta */
            _qadic_dense_artin_schreier_root_ii(Delta, Delta_y, Delta_x, V, roots, m - k, b + k, qctx + k);
            Delta->val += a[i + 1];

            /* rop */
            qadic_dense_add(rop, rop, Delta, qctxi);
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

/* Newton lift for modular polynomials */
/* k = 0, M = N' */
void qadic_dense_gen_newton_lift_ii_modular(qadic_dense_t rop, const qadic_dense_t op,
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

        _qadic_dense_gen_newton_lift_ii_modular(rop, op, roots, m, b, qctx);

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

    fmpz_t two, pow;
    padic_t ptwo, r, rinv;

    padic_ctx_t pctx2, pctx3;
    padic_poly_t modulus;
    qadic_dense_ctx_t ctx2, ctx3;

    qadic_dense_t gamma, gamma2, t, w, z;

    fmpz_init(two);
    fmpz_set_ui(two, 2);
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
    qadic_dense_zero(w);
    for (i = v; i < N + s - 1 - 2*v; i += 2*v)
    {
        padic_inv(rinv, r, pctx3);
        padic_poly_scalar_mul_padic(t, gamma, rinv, pctx3);
        padic_poly_reduce(t, pctx3);
        qadic_dense_add(w, w, t, ctx3);
        qadic_dense_mul(gamma, gamma, gamma2, ctx3);
        padic_add(r, r, ptwo, pctx3);
    }
    {
        padic_inv(rinv, r, pctx3);
        padic_poly_scalar_mul_padic(t, gamma, rinv, pctx3);
        padic_poly_reduce(t, pctx3);
        qadic_dense_add(w, w, t, ctx3);
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

    fmpz_clear(two);
    fmpz_clear(pow);
}

int main(int argc, char *argv[])
{
    clock_t time, total_time;
    long d, N;
    fmpz_t p, pow, s, t2;
    padic_t t;
    padic_struct *traces;
    fmpz_poly_t f, g;
    padic_poly_t h;
    padic_ctx_t pctx;
    qadic_dense_ctx_t qctx, qctx_one;
    qadic_dense_struct *roots;
    qadic_dense_t a, x, xinv;

    total_time = clock();

    if (argc < 2)
    {
        printf("Usage: points d [N]\n");
        return 0;
    }

    d = atol(argv[1]);

    if (argc > 2)
        N = atol(argv[2]);
    else
        N = ((d + 1) >> 1 ) + 2;

    fmpz_init(p);
    fmpz_set_ui(p, 2);

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
    time = clock();
    */
    qadic_dense_traces_init(&traces, qctx);
    printf("Traces took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);

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
    qadic_dense_gen_newton_lift_ii_modular(x, a, roots, qctx);
#if DEBUG >= 1
    printf("x = ");
    qadic_dense_print_pretty(x, qctx);
    printf("\n");
#endif
    printf("Lift took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);

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
    padic_init(t, pctx);
    qadic_dense_norm_char_2(t, xinv, traces, qctx);
    /* qadic_dense_norm(t, xinv, qctx); */
    printf("norm = ");
    padic_print(t, pctx);
    printf("\n");
    printf("Norm took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);

    fmpz_init(t2);
    padic_get_fmpz(t2, t, pctx);
    fmpz_init(s);
    fmpz_pow_ui(s, t2, 2);
    fmpz_init(pow);
    fmpz_pow_ui(pow, p, d + 1);

    if (fmpz_cmp(s, pow) > 0)
    {
        fmpz_pow_ui(pow, p, N);
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

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);

    fmpz_clear(p);
    fmpz_clear(pow);
    fmpz_clear(s);
    fmpz_clear(t2);

    return 0;
}
