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
#include "qadic.h"
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

#define MGM_LENGTH 23
long polymgm[MGM_LENGTH][9] = {
    {240, 8, 5, 3, 0, 0, 0, 0, 0},
    {444, 9, 7, 5, 4, 3, 2, 1, 0},
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
};

void fmpz_poly_init_ui(fmpz_poly_t rop, long* op)
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

void qadic_ctx_modulus(padic_poly_t rop, const qadic_ctx_t ctx)
{
    long i;

    padic_poly_fit_length(rop, ctx->len);
    padic_poly_zero(rop);

    for (i = 0; i < ctx->len; i++)
        fmpz_set(rop->coeffs + ctx->j[i], ctx->a + i);

    _padic_poly_set_length(rop, ctx->len);
    padic_poly_canonicalise(rop, (&ctx->pctx)->p);
}

void qadic_ctx_init(qadic_ctx_t ctx, const padic_poly_t modulus,
                    const padic_ctx_t pctx, const char *var)
{
    long d;
    long i, j;

    /* Find number of non-zero coefficients */
    d = padic_poly_degree(modulus);
    ctx->len = 1;

    for (i = 0; i < d; i++)
    {
        /* if (!padic_is_zero(padic_poly_get_coeff_padic(coeff, modulus, i, pctx),
                              pctx)) */
        if (!fmpz_is_zero(modulus->coeffs + i))
            ctx->len ++;
    }

    ctx->a = _fmpz_vec_init(ctx->len);
    ctx->j = flint_malloc(ctx->len * sizeof(long));

    /* Copy the polynomial */
    j = 0;

    for (i = 0; i < d; i++)
    {
        if (!fmpz_is_zero(modulus->coeffs + i))
        {
            fmpz_set(ctx->a + j, modulus->coeffs + i);
            ctx->j[j] = i;
            j++;
        }
    }

    fmpz_set(ctx->a + j, modulus-> coeffs + d);
    ctx->j[j] = d;

    /* Complete the initialisation of the context */
    padic_ctx_init(&ctx->pctx, pctx->p, pctx->N, pctx->mode);

    ctx->var = flint_malloc(strlen(var));
    strcpy(ctx->var, var);

    return;
}

void qadic_ctx_init_reduce(qadic_ctx_t rop, const qadic_ctx_t op, long N)
{
    long i;

    fmpz_t pow;
    fmpz_init(pow);
    fmpz_pow_ui(pow, (&op->pctx)->p, N);

    rop->len = op->len;

    rop->a = _fmpz_vec_init(rop->len);
    rop->j = flint_malloc(rop->len * sizeof(long));

    for (i = 0; i < rop->len; i++)
    {
        rop->j[i] = op->j[i];
        fmpz_mod(rop->a + i, op->a + i, pow);
    }

    padic_ctx_init(&rop->pctx, (&op->pctx)->p, N, (&op->pctx)->mode);

    rop->var = flint_malloc(strlen(op->var));
    strcpy(rop->var, op->var);

    fmpz_clear(pow);
}

int fmpz_poly_init_conway(fmpz_poly_t poly,
                          const fmpz_t p, long d)
{
    char *buf;
    FILE *file;

    if (fmpz_cmp_ui(p, 109987) > 0)
    {
        printf("Exception (fmpz_poly_init_conway).  Conway polynomials \n");
        printf("are only available for primes up to 109987.\n");
        return 1;
    }

    buf  = flint_malloc(832);
    file = fopen(FLINT_CPIMPORT, "r");

    if (!file)
    {
        printf("Exception (fmpz_poly_init_conway).  File loading.\n");
        return 2;
    }

    while (fgets(buf, 832, file))
    {
        char *tmp = buf;

        /* Different prime? */
        if (fmpz_cmp_ui(p, atoi(tmp)))
            continue;

        while (*tmp++ != ' ') ;

        /* Same degree? */
        if (d == atoi(tmp))
        {
            long i;
            char *ptr;

            /* Initialization */
            fmpz_poly_init2(poly, d + 1);

            /* Read coefficients */
            ptr = tmp;

            for (i = 0; i <= d; i++)
            {
                while (*ptr++ != ' ') ;

                fmpz_poly_set_coeff_si(poly, i, atoi(ptr));
            }

            fclose(file);
            flint_free(buf);
            return 0;
        }
    }

    fclose(file);
    flint_free(buf);

    printf("Exception (fmpz_poly_init_conway).  The polynomial for \n");
    printf("(p,d) = (%ld,%ld) is not present in the database.\n", *p, d);
    return 1;
}

int padic_poly_init_conway(padic_poly_t poly,
                           const fmpz_t p, long d,
                           const padic_ctx_t ctx)
{
    fmpz_poly_t cpoly;
    if (fmpz_poly_init_conway(cpoly, p, d))
    {
        return 1;
    }
    padic_poly_init2(poly, d + 1);
    padic_poly_set_fmpz_poly(poly, cpoly, ctx);
    fmpz_poly_clear(cpoly);
    return 0;
}

/* Powers of p should be shared */
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
        long *a, i, j, l, l0, l1, ld, n, t;
        fmpz *pow;
        fmpz_poly_t delta0, df0, delta1, df1, W, Delta;

        /* Degree and length */
        l = d + 1;
        l0 = (l + 1) >> 1;
        l1 = (l >> 1);

        /* Number of steps */
        n = FLINT_CLOG2(N) + 1;

        /* Precision at each step */
        a = flint_malloc(n * sizeof(long));

        for (a[i = 0] = N; a[i] > 1; i++)
            a[i + 1] = (a[i] + 1) / 2;

        pow = _fmpz_vec_init(n);

        fmpz_poly_init2(delta0, l0);
        fmpz_poly_init2(df0, l);
        fmpz_poly_init2(delta1, l1);
        fmpz_poly_init2(df1, l);
        fmpz_poly_init2(W, l);
        fmpz_poly_init2(Delta, l);

        /* Compute powers of 2 */
        {
            t = 0;
            fmpz_set(pow + i, p);
        }
        for (i--; i >= 1; i--)
        {
            if (a[i] & 1L)
            {
                fmpz_mul_2exp(pow + i, pow + (i + 1), t);
                t <<= 1;
            }
            else
            {
                t += a[i + 1];
                fmpz_mul_2exp(pow + i, pow + (i + 1), a[i + 1]);
            }
        }
        {
            if (a[i] & 1L)
                fmpz_mul_2exp(pow + i, pow + (i + 1), t);
            else
                fmpz_mul_2exp(pow + i, pow + (i + 1), a[i + 1]);
        }

        /* Lifting */
        i = n - 1;
        {
            pow = pow + i;

            fmpz_poly_set(delta, V);

            ld = delta->length;

            _fmpz_vec_scalar_mod_fmpz(delta->coeffs, delta->coeffs, ld, pow);
            _fmpz_poly_normalise(delta);
        }
        for (i--; i >= 0; i--)
        {
            pow--;

            ld = delta->length;

            /* delta_0, delta_1, and products */
            _fmpz_poly_set_length(delta0, (ld + 1) >> 1);
            _fmpz_poly_set_length(delta1, (ld >> 1));

            for (j = 0; j < ld; j++)
                fmpz_set(((j % 2)?delta1:delta0)->coeffs + (j >> 1), delta->coeffs + j);

            _fmpz_poly_normalise(delta0);
            _fmpz_poly_normalise(delta1);

            _fmpz_vec_scalar_mod_fmpz(df0->coeffs, f0->coeffs, f0->length, pow);
            _fmpz_vec_scalar_mod_fmpz(df1->coeffs, f1->coeffs, f1->length, pow);
            _fmpz_poly_set_length(df0, f0->length);
            _fmpz_poly_set_length(df1, f1->length);
            _fmpz_poly_normalise(df0);
            _fmpz_poly_normalise(df1);

            _fmpz_mod_poly_mul(df0->coeffs, delta0->coeffs, delta0->length, f0->coeffs, f0->length, pow);
            _fmpz_mod_poly_mul(df1->coeffs, delta1->coeffs, delta1->length, f1->coeffs, f1->length, pow);

            _fmpz_poly_set_length(df0, delta0->length + f0->length - 1);
            _fmpz_poly_set_length(df1, delta1->length + f1->length);

            _fmpz_poly_shift_left(df1->coeffs, df1->coeffs, df1->length - 1, 1);

            _fmpz_mod_poly_sub(df0->coeffs, df0->coeffs, df0->length, df1->coeffs, df1->length, pow);

            _fmpz_poly_set_length(df0, FLINT_MAX(df0->length, df1->length));
            _fmpz_poly_normalise(df0);

            /* Use vec method */
            for (j = 0; j < df0->length; j++)
                fmpz_mul_2exp(df0->coeffs + j, df0->coeffs + j, 1);
            /* Computed everything above with one superfluous bit... */
            _fmpz_vec_scalar_mod_fmpz(df0->coeffs, df0->coeffs, df0->length, pow);
            _fmpz_poly_normalise(df0);

            /* W */
            _fmpz_mod_poly_add(W->coeffs, delta->coeffs, ld, V->coeffs, V->length, pow);
            _fmpz_poly_set_length(W, FLINT_MAX(ld, V->length));
            _fmpz_poly_normalise(W);

            _fmpz_mod_poly_sub(W->coeffs, W->coeffs, W->length, df0->coeffs, df0->length, pow);
            _fmpz_poly_set_length(W, FLINT_MAX(W->length, df0->length));
            _fmpz_poly_normalise(W);

            /* Use vec method */
            for (j = 0; j < W->length; j++)
                fmpz_fdiv_q_2exp(W->coeffs + j, W->coeffs + j, a[i + 1]);

            /* Delta */
            _padic_poly_teichmuller_inc(Delta, f0, f1, W, d, p, a[i] - a[i + 1]);

            /* Use vec method */
            for (j = 0; j < Delta->length; j++)
                fmpz_mul_2exp(Delta->coeffs + j, Delta->coeffs + j, a[i + 1]);
            /* Check precision here */
            _fmpz_vec_scalar_mod_fmpz(Delta->coeffs, Delta->coeffs, Delta->length, pow);
            _fmpz_poly_normalise(Delta);

            /* update delta */
            _fmpz_mod_poly_add(delta->coeffs, delta->coeffs, ld, Delta->coeffs, Delta->length, pow);
            _fmpz_poly_set_length(delta, FLINT_MAX(ld, Delta->length));
            _fmpz_poly_normalise(delta);
        }

        fmpz_poly_clear(delta0);
        fmpz_poly_clear(df0);
        fmpz_poly_clear(delta1);
        fmpz_poly_clear(df1);
        fmpz_poly_clear(W);
        fmpz_poly_clear(Delta);

        _fmpz_vec_clear(pow, n);

        flint_free(a);
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
    else if (*(p) == 2L)
    {
        long *a, d, i, j, l, l0, l1, n, t;
        /* One should be sufficient */
        fmpz *pow;
        fmpz_poly_t f0, f0_sqr, f1, f1_sqr, V, delta;

        /* Degree and length */
        d = fmpz_poly_degree(f);
        l = d + 1;
        l0 = (l + 1) >> 1;
        l1 = (l >> 1);

        /* Number of steps */
        n = FLINT_CLOG2(N) + 1;

        /* Precision at each step */
        a = flint_malloc(n * sizeof(long));

        for (a[i = 0] = N; a[i] > 1; i++)
            a[i + 1] = (a[i] + 1) / 2;

        pow = _fmpz_vec_init(n);

        fmpz_poly_init2(f0, l0);
        fmpz_poly_init2(f0_sqr, l);
        fmpz_poly_init2(f1, l1);
        fmpz_poly_init2(f1_sqr, l);
        fmpz_poly_init2(V, l);
        fmpz_poly_init2(delta, l);

        /* Compute powers of 2 */
        {
            t = 0;
            fmpz_set(pow + i, p);
        }
        for (i--; i >= 1; i--)
        {
            if (a[i] & 1L)
            {
                fmpz_mul_2exp(pow + i, pow + (i + 1), t);
                t <<= 1;
            }
            else
            {
                t += a[i + 1];
                fmpz_mul_2exp(pow + i, pow + (i + 1), a[i + 1]);
            }
        }
        {
            if (a[i] & 1L)
                fmpz_mul_2exp(pow + i, pow + (i + 1), t);
            else
                fmpz_mul_2exp(pow + i, pow + (i + 1), a[i + 1]);
        }

        /* Lifting */
        i = n - 1;
        {
            fmpz_poly_set(g, f);

            pow = pow + i;
        }
        for (i--; i >= 0; i--)
        {
            pow--;

            /* f_0, f_1, and squares */
            _fmpz_poly_set_length(f0, l0);
            _fmpz_poly_set_length(f1, l1);

            for (j = 0; j < l; j++)
                fmpz_set(((j % 2)?f1:f0)->coeffs + (j >> 1), g->coeffs + j);

            _fmpz_poly_normalise(f0);
            _fmpz_poly_normalise(f1);

            _fmpz_mod_poly_sqr(f0_sqr->coeffs, f0->coeffs, f0->length, pow);
            _fmpz_mod_poly_sqr(f1_sqr->coeffs, f1->coeffs, f1->length, pow);

            _fmpz_poly_set_length(f0_sqr, (f0->length << 1) - 1);
            _fmpz_poly_set_length(f1_sqr, (f1->length << 1));

            _fmpz_poly_shift_left(f1_sqr->coeffs, f1_sqr->coeffs, f1_sqr->length - 1, 1);

            /* V */
            _fmpz_poly_set_length(V, l);

            /* Use vec method */
            for (j = 0; j < l; j++)
                fmpz_set(V->coeffs + j, g->coeffs + j);

            _fmpz_mod_poly_sub(V->coeffs, V->coeffs, V->length, f0_sqr->coeffs, f0_sqr->length, pow);
            _fmpz_poly_normalise(V);

            _fmpz_mod_poly_add(V->coeffs, V->coeffs, V->length, f1_sqr->coeffs, f1_sqr->length, pow);
            _fmpz_poly_set_length(V, FLINT_MAX(V->length, f1_sqr->length));
            _fmpz_poly_normalise(V);

            /* Use vec method */
            for (j = 0; j < V->length; j++)
                fmpz_fdiv_q_2exp(V->coeffs + j, V->coeffs + j, a[i + 1]);

            /* delta */
            _padic_poly_teichmuller_inc(delta, f0, f1, V, d, p, a[i] - a[i + 1]);

            /* Use vec method */
            for (j = 0; j < delta->length; j++)
                fmpz_mul_2exp(delta->coeffs + j, delta->coeffs + j, a[i + 1]);
            _fmpz_vec_scalar_mod_fmpz(delta->coeffs, delta->coeffs, delta->length, pow);
            _fmpz_poly_normalise(delta);

            /* update g */
            _fmpz_mod_poly_add(g->coeffs, g->coeffs, l, delta->coeffs, delta->length, pow);
        }

        if (!fmpz_is_one(fmpz_poly_lead(g)))
            _fmpz_mod_poly_neg(g->coeffs, g->coeffs, g->length, pow);

        fmpz_poly_clear(f0);
        fmpz_poly_clear(f0_sqr);
        fmpz_poly_clear(f1);
        fmpz_poly_clear(f1_sqr);
        fmpz_poly_clear(V);
        fmpz_poly_clear(delta);

        _fmpz_vec_clear(pow, n);

        flint_free(a);
    }
}

/* This really takes place in a finite field of small characteristic */
void _qadic_proot_basis_init(qadic_struct *roots, const qadic_ctx_t ctx)
{
    long j, p, d;
    qadic_t t, r;
    fmpz_t pow;
    d = qadic_ctx_degree(ctx);
    p = *(&ctx->pctx)->p;

    fmpz_init(pow);
    fmpz_pow_ui(pow, &p, d - 1);

    padic_poly_init2(t, d);
    _padic_poly_set_length(t, 2);
    fmpz_set_ui(t->coeffs, 0);
    fmpz_set_ui(t->coeffs + 1, 1);
    qadic_pow(t, t, pow, ctx);

    qadic_init(r);
    qadic_one(r, ctx);

    j = 0;
    {
        qadic_set(roots, r);
    }
    for (j++; j < p; j++)
    {
        qadic_mul(r, r, t, ctx);
        qadic_set(roots + j, r);
    }

    qadic_clear(t);
    qadic_clear(r);
}

void qadic_proot_basis_init(qadic_struct **roots, const qadic_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;
    const long p = *(&ctx->pctx)->p;

    if (N != 1)
    {
        printf("Error (qadic_init_proot_basis).  precision different from one\n");
        abort();
    }
    else if (COEFF_IS_MPZ(p))
    {
        printf("Error (qadic_init_proot_basis).  characteristic does not fit in a long\n");
        abort();
    }
    else
    {
        long i;
        *roots = flint_malloc(p * sizeof(qadic_struct));
        for (i = 0; i < p; i++)
            qadic_init(*roots + i);
        _qadic_proot_basis_init(*roots, ctx);
    }
}

void _qadic_proot(qadic_t rop, const qadic_t op, const qadic_struct *roots,
                  const qadic_ctx_t ctx)
{
    const long d = qadic_ctx_degree(ctx);
    const long p = *(&ctx->pctx)->p;
    long j, k;
    qadic_t t;

    qadic_init(t);

    qadic_zero(rop);

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

        qadic_mul(t, t, roots + j, ctx);

        qadic_add(rop, rop, t, ctx);
    }

    qadic_clear(t);
}

void qadic_proot(qadic_t rop, const qadic_t op, const qadic_struct *roots,
                 const qadic_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;
    const long p = *(&ctx->pctx)->p;

    if (N != 1)
    {
        printf("ERROR (qadic_proot).  not implemented.\n");
        abort();
    }
    else if (COEFF_IS_MPZ(p))
    {
        printf("Error (qadic_proot).  characteristic does not fit in a long\n");
        abort();
    }
    else if (qadic_is_zero(op))
    {
        qadic_zero(rop);
        return;
    }
    else
    {
        if (rop == op)
        {
            qadic_t t;
            qadic_init(t);
            qadic_set(t, op);
            _qadic_proot(rop, t, roots, ctx);
            qadic_clear(t);
        }
        else
        {
            _qadic_proot(rop, op, roots, ctx);
        }
    }
}

void qadic_frobenius_teichmuller_init(qadic_struct **frobs,
                                      const qadic_ctx_t ctx)
{
    long i, d;
    qadic_t t;

    d = qadic_ctx_degree(ctx);

    padic_poly_init2(t, d - 1);
    _padic_poly_set_length(t, 2);
    fmpz_set_ui(t->coeffs, 0);
    fmpz_set_ui(t->coeffs + 1, 1);
    qadic_pow(t, t, (&ctx->pctx)->p, ctx);

    *frobs = flint_malloc(d * sizeof(qadic_struct));

    i = 0;
    {
        qadic_init(*frobs);
        qadic_one(*frobs, ctx);
    }
    for (i++; i < d; i++)
    {
        qadic_init(*frobs + i);
        qadic_mul(*frobs + i, *frobs + i - 1, t, ctx);
    }
}


void _qadic_frobenius_teichmuller(qadic_t rop, const qadic_t op,
                                 const qadic_struct *frobs,
                                 const qadic_ctx_t ctx)
{
    long i, d;
    qadic_t t;
    padic_t c;

    d = qadic_ctx_degree(ctx);

    padic_init(c, &ctx->pctx);
    qadic_init(t);

    qadic_zero(rop);

    for (i = 0; i < d; i++)
    {
        padic_poly_set(t, frobs + i);
        padic_poly_reduce(t, &ctx->pctx);

        padic_poly_get_coeff_padic(c, op, i, &ctx->pctx);
        padic_poly_scalar_mul_padic(t, t, c, &ctx->pctx);

        qadic_add(rop, rop, t, ctx);
    }
}

void qadic_frobenius_teichmuller(qadic_t rop, const qadic_t op,
                                 const qadic_struct *frobs,
                                 const qadic_ctx_t ctx)
{
    if (rop == op)
    {
        qadic_t t;
        qadic_set(t, op);
        _qadic_frobenius_teichmuller(rop, t, frobs, ctx);
        qadic_clear(t);
    }
    else
    {
        _qadic_frobenius_teichmuller(rop, op, frobs, ctx);
    }
}

void qadic_artin_schreier_root_ii(qadic_t x, const qadic_t alpha,
                                   const qadic_t beta, const qadic_t gamma,
                                   const qadic_struct *roots,
                                   const qadic_struct *frobs,
                                   const qadic_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (N == 1)
    {
        qadic_inv(x, alpha, ctx);
        qadic_mul(x, x, gamma, ctx);
        qadic_neg(x, x, ctx);
        qadic_proot(x, x, roots, ctx);
    }
    else
    {
        long *a, i, n, M;
        qadic_t alpha2, beta2, gamma2, gamma3, Delta2;
        qadic_ctx_t ctx2;
        qadic_ctx_struct *qctx;

        /* Number of steps */
        n = FLINT_CLOG2(N) + 1;

        /* Precision at each step */
        a = flint_malloc(n * sizeof(long));

        for (a[i = 0] = N; a[i] > 1; i++)
            a[i + 1] = (a[i] + 1) / 2;

        /* Contexts */
        qctx = flint_malloc(n * sizeof(qadic_ctx_struct));

        for (i = 0; i < n; i++)
            qadic_ctx_init_reduce(qctx + i, ctx, a[i]);

        qadic_init(alpha2);
        qadic_init(beta2);
        qadic_init(gamma2);
        qadic_init(gamma3);
        qadic_init(Delta2);

        /* Lifting */
        i = n - 1;
        {
            qctx = qctx + i;

            qadic_set(alpha2, alpha);
            qadic_set(beta2, beta);
            qadic_set(gamma2, gamma);

            padic_poly_reduce(alpha2, &qctx->pctx);
            padic_poly_reduce(beta2, &qctx->pctx);
            padic_poly_reduce(gamma2, &qctx->pctx);

            qadic_inv(x, alpha2, qctx);
            qadic_mul(x, x, gamma2, qctx);
            qadic_neg(x, x, qctx);
            qadic_proot(x, x, roots, qctx);
        }
        for (i--; i >= 0; i--)
        {
            qctx--;
            M = a[i] - a[i + 1];

            qadic_set(alpha2, alpha);
            qadic_set(beta2, beta);
            qadic_set(gamma2, gamma);

            padic_poly_reduce(alpha2, &qctx->pctx);
            padic_poly_reduce(beta2, &qctx->pctx);
            padic_poly_reduce(gamma2, &qctx->pctx);

            qadic_frobenius_teichmuller(gamma3, x, frobs, qctx);
            qadic_mul(gamma3, gamma3, alpha2, qctx);
            qadic_mul(Delta2, beta2, x, qctx);
            qadic_add(gamma3, gamma3, Delta2, qctx);
            qadic_add(gamma3, gamma3, gamma2, qctx);
            gamma3->val -= a[i + 1];

            qadic_ctx_init_reduce(ctx2, qctx, M);

            padic_poly_reduce(alpha2, &ctx2->pctx);
            padic_poly_reduce(beta2, &ctx2->pctx);
            padic_poly_reduce(gamma3, &ctx2->pctx);

            qadic_artin_schreier_root_ii(Delta2, alpha2, beta2, gamma3, roots, frobs, ctx2);
            Delta2->val += a[i + 1];

            qadic_add(x, x, Delta2, qctx);

            qadic_ctx_clear(ctx2);
        }

        qadic_clear(alpha2);
        qadic_clear(beta2);
        qadic_clear(gamma2);
        qadic_clear(gamma3);
        qadic_clear(Delta2);

        for (i = 0; i < n; i++)
            qadic_ctx_clear(qctx + i);

        flint_free(qctx);

        flint_free(a);
    }
}

/* Newton lift for modular polynomials */
/* k = 0, M = N' */
void qadic_gen_newton_lift_ii_modular(qadic_t rop, const qadic_t op,
                                      const qadic_struct *roots,
                                      const qadic_struct *frobs,
                                      const qadic_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (N == 1)
    {
        qadic_set(rop, op);
    }
    else if (*(&ctx->pctx)->p == 2L)
    {
        long *a, i, n;
        qadic_t Delta_x, Delta_y, Delta, V, y, xy;
        qadic_ctx_struct *qctx;

        qadic_init(Delta_x);
        qadic_init(Delta_y);
        qadic_init(Delta);
        qadic_init(V);
        qadic_init(y);
        qadic_init(xy);

        /* Number of steps */
        n = FLINT_CLOG2(N) + 1;

        /* Precision at each step */
        a = flint_malloc(n * sizeof(long));

        for (a[i = 0] = N; a[i] > 1; i++)
            a[i + 1] = (a[i] + 1) / 2;

        /* Contexts */
        qctx = flint_malloc(n * sizeof(qadic_ctx_struct));

        for (i = 0; i < n; i++)
            qadic_ctx_init_reduce(qctx + i, ctx, a[i]);

        /* Lifting */
        i = n - 1;
        {
            qctx = qctx + i;

            qadic_set(rop, op);
        }
        for (i--; i >= 0; i--)
        {
            qctx--;

            /* y */
            qadic_frobenius_teichmuller(y, rop, frobs, qctx);
            /* 2 y */
            qadic_set(Delta_x, y);
            Delta_x->val += 1;
            padic_poly_reduce(Delta_x, &qctx->pctx);
            /* 2 x y */
            qadic_mul(xy, rop, Delta_x, qctx);
            /* x + 2 y */
            qadic_add(Delta, rop, Delta_x, qctx);
            /* 8 x y */
            qadic_set(Delta_x, xy);
            Delta_x->val += 2;
            padic_poly_reduce(Delta_x, &qctx->pctx);
            /* x + 2 y + 8 x y */
            qadic_add(Delta, Delta, Delta_x, qctx);

            /* V */
            /* (x + 2 y + 8 x y)² */
            qadic_pow(V, Delta, (&qctx->pctx)->p, qctx);
            /* (x + 2 y + 8 x y)² + y */
            qadic_add(V, V, y, qctx);
            /* 4 x y */
            xy->val += 1;
            padic_poly_reduce(xy, &qctx->pctx);
            /* (x + 2 y + 8 x y)² + y + 4 x y */
            qadic_add(V, V, xy, qctx);
            /* (x + 2 y + 8 x y)² + y + 4 x y / 2^N' */
            V->val -= a[i + 1];

            /* Delta_x */
            /* x + 2 y + 8 x y */
            padic_poly_reduce(Delta, &(qctx + 1)->pctx);
            /* 1 */
            qadic_one(xy, qctx + 1);
            /* 8 y */
            qadic_set(Delta_x, y);
            Delta_x->val += 3;
            padic_poly_reduce(Delta_x, &(qctx + 1)->pctx);
            /* 1 + 8 y */
            qadic_add(xy, xy, Delta_x, qctx + 1);
            /* (x + 2 y + 8 x y) (1 + 8 y) */
            qadic_mul(Delta_x, Delta, xy, qctx + 1);
            /* 2 (x + 2 y + 8 x y) (1 + 8 y) */
            Delta_x->val += 1;
            padic_poly_reduce(Delta_x, &(qctx + 1)->pctx);
            /* 4 y */
            y->val += 2;
            padic_poly_reduce(y, &(qctx + 1)->pctx);
            /* 2 (x + 2 y + 8 x y) (1 + 8 y) + 4 y*/
            qadic_add(Delta_x, Delta_x, y, qctx + 1);

            /* Delta_y */
            /* (x + 2 y + 8 x y) */
            Delta->val += 2;
            padic_poly_reduce(Delta, &(qctx + 1)->pctx);
            /* 1 + 2 (x + 2 y + 8 x y) */
            qadic_one(xy, qctx + 1);
            qadic_add(Delta_y, Delta, xy, qctx + 1);
            /* 4 x */
            qadic_set(y, rop);
            y->val += 2;
            padic_poly_reduce(y, &(qctx + 1)->pctx);
            /* 1 + 4 x */
            qadic_add(xy, xy, y, qctx + 1);
            /* (1 + 2 (x + 2 y + 8 x y)) (1 + 4 x) */
            qadic_mul(Delta_y, Delta_y, xy, qctx + 1);

            /* Delta */
            qadic_artin_schreier_root_ii(Delta, Delta_y, Delta_x, V, roots, frobs, qctx + 1);
            Delta->val += a[i + 1];

            /* rop */
            qadic_add(rop, rop, Delta, qctx);
        }

        qadic_clear(Delta_x);
        qadic_clear(Delta_y);
        qadic_clear(Delta);
        qadic_clear(V);
        qadic_clear(y);
        qadic_clear(xy);

        for (i = 0; i < n; i++)
            qadic_ctx_clear(qctx + i);

        flint_free(qctx);

        flint_free(a);
    }
    else
    {
        printf("ERROR (_padic_poly_teichmuller).  odd characteristic not implemented.\n");
        abort();
    }
}

void qadic_norm_char_2(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;
    const long s = floor(sqrt(N)) / 2;

    long i, v;

    fmpz_t two, pow;
    padic_t ptwo, r, rinv;

    padic_ctx_t pctx2, pctx3;
    padic_poly_t modulus;
    qadic_ctx_t ctx2, ctx3;

    qadic_t gamma, gamma2, t, w, z;

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
    qadic_ctx_modulus(modulus, ctx);
    qadic_ctx_init(ctx2, modulus, pctx2, ctx->var);
    qadic_ctx_init(ctx3, modulus, pctx3, ctx->var);

    qadic_init(gamma);
    qadic_init(gamma2);
    qadic_init(t);
    qadic_init(w);
    qadic_init(z);

    /* z = (a^(2^s) - 1) / 2 */
    qadic_pow(z, op, pow, ctx2);
    qadic_one(w, ctx2);
    qadic_sub(z, z, w, ctx2);
    z->val -= 1;
    v = z->val;

    /* gamma = z / (1 + z)*/
    qadic_add(gamma, z, w, ctx3);
    qadic_inv(gamma, gamma, ctx3);
    qadic_mul(gamma, gamma, z, ctx3);
    /* gamma2 = gamma^2*/
    qadic_pow(gamma2, gamma, two, ctx3);

    /* w = 2^(-s) log(a) */
    padic_one(r, pctx3);
    qadic_zero(w);
    for (i = v; i < N + s - 1 - 2*v; i += 2*v)
    {
        padic_inv(rinv, r, pctx3);
        padic_poly_scalar_mul_padic(t, gamma, rinv, pctx3);
        padic_poly_reduce(t, pctx3);
        qadic_add(w, w, t, ctx3);
        qadic_mul(gamma, gamma, gamma2, ctx3);
        padic_add(r, r, ptwo, pctx3);
    }
    {
        padic_inv(rinv, r, pctx3);
        padic_poly_scalar_mul_padic(t, gamma, rinv, pctx3);
        padic_poly_reduce(t, pctx3);
        qadic_add(w, w, t, ctx3);
    }
    w->val -= s - 1;

    /* tr(w) */
    /* use precomputations? */
    qadic_trace(rop, w, ctx);

    /* exp(tr(w)) */
    padic_exp(rop, rop, &ctx->pctx);

    qadic_clear(gamma);
    qadic_clear(gamma2);
    qadic_clear(t);
    qadic_clear(w);
    qadic_clear(z);

    padic_ctx_clear(pctx2);
    padic_ctx_clear(pctx3);
    padic_poly_clear(modulus);
    qadic_ctx_clear(ctx2);
    qadic_ctx_clear(ctx3);

    padic_clear(ptwo, &ctx->pctx);
    padic_clear(r, &ctx->pctx);
    padic_clear(rinv, &ctx->pctx);

    fmpz_clear(two);
    fmpz_clear(pow);
}

int main(int argc, char* argv[])
{
    clock_t time, total_time;
    long d, N;
    fmpz_t p, pow, s, t2;
    fmpz_poly_t f, g;
    padic_poly_t h;
    padic_ctx_t pctx;
    qadic_ctx_t qctx, qctx_one;
    qadic_struct *frobs, *roots;
    qadic_t a, x, xinv;
    padic_t t;

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


    qadic_ctx_init(qctx, h, pctx, "t");
    qadic_ctx_init_reduce(qctx_one, qctx, 1);

    time = clock();
    qadic_proot_basis_init(&roots, qctx_one);
    printf("p-th root took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
    time = clock();
    qadic_frobenius_teichmuller_init(&frobs, qctx);
    printf("Frobenius took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);
#if DEBUG >= 1
    printf("roots = ");
    qadic_print_pretty(roots + 1, qctx);
    printf("\n");
    printf("frobs = ");
    qadic_print_pretty(frobs + 1, qctx);
    printf("\n");
#endif

    padic_poly_init2(a, d - 1);
    _padic_poly_set_length(a, d - 1);
    fmpz_set_ui(a->coeffs + d - 2, 1);
    printf("a = ");
    qadic_print_pretty(a, qctx);
    printf("\n");

    time = clock();
    qadic_init(x);
    qadic_gen_newton_lift_ii_modular(x, a, roots, frobs, qctx);
#if DEBUG >= 1
    printf("x = ");
    qadic_print_pretty(x, qctx);
    printf("\n");
#endif
    printf("Lift took %f seconds.\n", ((double) (clock() - time)) / CLOCKS_PER_SEC);

    time = clock();
    qadic_init(xinv);
    qadic_set(xinv, x);
    /* groumpf, that means, we can gain 2 bits above... */
    xinv->val += 2;
    padic_poly_reduce(xinv, &qctx->pctx);
    qadic_one(x, qctx);
    qadic_add(xinv, xinv, x, qctx);
    qadic_inv(xinv, xinv, qctx);
#if DEBUG >= 1
    printf("xinv = ");
    qadic_print_pretty(xinv, qctx);
    printf("\n");
#endif

    /* Go back to sparse modulus from here... */
    padic_init(t, pctx);
    qadic_norm_char_2(t, xinv, qctx);
    /* qadic_norm(t, xinv, qctx); */
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

    qadic_clear(a);
    qadic_clear(x);
    qadic_clear(xinv);

    qadic_ctx_clear(qctx);
    padic_ctx_clear(pctx);

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);

    fmpz_clear(p);
    fmpz_clear(pow);
    fmpz_clear(s);
    fmpz_clear(t2);

    return 0;
}
