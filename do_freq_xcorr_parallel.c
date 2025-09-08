/*  $Id: do_freq_xcorr_parallel.c 1 2025-08-15Z parallel $  */
/*  Parallel version: Frequency domain cross-correlation  */

#include "gmt.h"
#include "gmtsar.h"
#include "xcorr_parallel.h"
#include <math.h>
#include <string.h>
#include <omp.h>

/* Thread-safe FFT multiplication (each thread uses independent buffers) */
static void fft_multiply_parallel(void *API, int N, int M,
                                  struct FCOMPLEX *c1,
                                  struct FCOMPLEX *c2,
                                  struct FCOMPLEX *c3) {
    int i, j, isign;

    /* Forward FFT */
    GMT_FFT_2D(API, (float *)c1, N, M, GMT_FFT_FWD, GMT_FFT_COMPLEX);
    GMT_FFT_2D(API, (float *)c2, N, M, GMT_FFT_FWD, GMT_FFT_COMPLEX);

    /* Multiply c1 with conj(c2) */
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            isign = (int)lrint(pow(-1.0, (double)(i + j)));
            c3[i * N + j] = RCmul(isign,
                                  Cmul(c1[i * N + j],
                                       Conjg(c2[i * N + j])));
        }
    }

    /* Inverse FFT for correlation matrix */
    GMT_FFT_2D(API, (float *)c3, N, M, GMT_FFT_INV, GMT_FFT_COMPLEX);
}

/* Main function: Frequency domain cross-correlation (parallel version) */
void do_freq_corr_parallel(void *API, struct xcorrp *xc, int iloc) {
    int i, j, ii;
    float ipeak = -999.0f, jpeak = -999.0f;
    double cmax = -1.0, cave = 0.0;
    double max_corr = -999.0;

    /* Get thread ID and local buffers */
    int tid = omp_get_thread_num();
    struct FCOMPLEX *c1 = xc->thread_data[tid].c1;
    struct FCOMPLEX *c2 = xc->thread_data[tid].c2;
    struct FCOMPLEX *c3 = xc->thread_data[tid].c3;

    /* Copy master data to thread buffers */
    memcpy(c1, xc->c1, xc->npx * xc->npy * sizeof(struct FCOMPLEX));
    memcpy(c2, xc->c2, xc->npx * xc->npy * sizeof(struct FCOMPLEX));

    /* FFT multiplication */
    fft_multiply_parallel(API, xc->npx, xc->npy, c1, c2, c3);

    /* Search for maximum correlation value and calculate average */
    #pragma omp parallel for collapse(2) reduction(+:cave) reduction(max:cmax)
    for (i = 0; i < xc->nyc; i++) {
        for (j = 0; j < xc->nxc; j++) {
            ii = (i + xc->ny_corr / 2) * xc->npx + j + (xc->nx_corr / 2);
            double val = Cabs(c3[ii]);
            cave += val;
            if (val > cmax) {
                cmax = val;
                /* Save indices */
                #pragma omp critical
                {
                    jpeak = j - xc->nxc / 2.0f;
                    ipeak = i - xc->nyc / 2.0f;
                }
            }
        }
    }

    if ((ipeak == -999.0f) || (jpeak == -999.0f)) {
        fprintf(stderr, "[PARALLEL] error! jpeak %f ipeak %f cmax %lf\n",
                jpeak, ipeak, cmax);
        exit(1);
    }

    /* Calculate average */
    cave /= (xc->nxc * xc->nyc);

    /* Initial maximum correlation (frequency domain estimation) */
    max_corr = (cmax / cave);

    /* Replace frequency domain estimation with time domain normalized correlation */
    if (ipeak != -999.0f) {
        max_corr = calc_time_corr_parallel(
            xc, (int)ipeak, (int)jpeak, xc->i1, xc->i2
        );
    }

    /* Scale correlation matrix */
    #pragma omp parallel for
    for (i = 0; i < xc->nxc * xc->nyc; i++) {
        int row = i / xc->nxc;
        int col = i % xc->nxc;
        int idx = (row + xc->ny_corr / 2) * xc->npx +
                  (col + xc->nx_corr / 2);
        xc->corr[i] = max_corr * Cabs(c3[idx]) / cmax;
    }

    /* Write results */
    xc->loc[iloc].xoff = -1 * jpeak;
    xc->loc[iloc].yoff = -1 * ipeak;
    xc->loc[iloc].corr = (float)max_corr;

    if (debug) {
        fprintf(stderr, " (freq) jpeak %f xoffset %d corr %4.2lf\n",
                jpeak, xc->x_offset, max_corr);
        fprintf(stderr, " (freq) ipeak %f yoffset %d corr %4.2lf\n",
                ipeak, xc->y_offset, max_corr);
    }
}