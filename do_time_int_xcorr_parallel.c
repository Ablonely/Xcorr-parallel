/*	$Id: do_time_int_xcorr_parallel.c 1 2023-10-01 00:00:00Z parallel $	*/
/*-------------------------------------------------------*/
#include "gmtsar.h"
#include "xcorr_parallel.h"
#include <math.h>
#include <omp.h>

/*-------------------------------------------------------------------------------*/
double calc_time_corr_parallel(struct xcorrp *xc, int ioff, int joff, int *i1, int *i2) {
    long long gamma_num = 0, gamma_denom1 = 0, gamma_denom2 = 0;
    double gamma, gamma_denom;
    const int npx = xc->npx;
    const int ysearch = xc->ysearch;
    const int xsearch = xc->xsearch;
    const int ny_corr = xc->ny_corr;
    const int nx_corr = xc->nx_corr;

    /* calculate normalized correlation */
    #pragma omp simd collapse(2) reduction(+:gamma_num, gamma_denom1, gamma_denom2)
    for (int ip = 0; ip < ny_corr; ip++) {
        for (int jp = 0; jp < nx_corr; jp++) {
            /* pixel values */
            const int idx1 = (ysearch + ip + ioff) * npx + jp + joff + xsearch;
            const int idx2 = (ysearch + ip) * npx + jp + xsearch;
            
            const long long a = i1[idx1];
            const long long b = i2[idx2];
            
            /* standard correlation */
            gamma_num += (a * b);
            gamma_denom1 += a * a;
            gamma_denom2 += b * b;
        }
    }

    gamma_denom = sqrt(1.0 * gamma_denom1 * gamma_denom2);

    if (gamma_denom == 0.0) {
        if (verbose)
            fprintf(stderr, "calc_corr: denominator = zero: setting corr to 0 \n");
        gamma = 0.0;
    } else {
        gamma = 100.0 * fabs(gamma_num / gamma_denom);
    }

    if (debug)
        fprintf(stdout, " corr %6.2lf \n", gamma);

    return gamma;
}

/*-------------------------------------------------------------------------------*/
double calc_time_corr_hat_parallel(struct xcorrp *xc, int ioff, int joff, int *i1, int *i2) {
    long long gamma_num = 0, gamma_denom1 = 0, gamma_denom2 = 0;
    double gamma_denom;
    double gamma;
    const int npx = xc->npx;
    const int ysearch = xc->ysearch;
    const int xsearch = xc->xsearch;
    const int ny_corr = xc->ny_corr;
    const int nx_corr = xc->nx_corr;

    #pragma omp simd collapse(2) reduction(+:gamma_num, gamma_denom1, gamma_denom2)
    for (int ip = 0; ip < ny_corr; ip++) {
        for (int jp = 0; jp < nx_corr; jp++) {
            /* pixel values */
            const int idx1 = (ysearch + ip + ioff) * npx + jp + joff + xsearch;
            const int idx2 = (ysearch + ip) * npx + jp + xsearch;
            
            const long long a = i1[idx1];
            const long long b = i2[idx2];
            const long long a2 = a * a;
            const long long b2 = b * b;
            
            /* frequency independent */
            gamma_num += (a2 * b2);
            gamma_denom1 += (a2 * a2);
            gamma_denom2 += (b2 * b2);
        }
    }

    gamma_denom = sqrtl(1.0 * gamma_denom1 * gamma_denom2);

    if (gamma_denom == 0.0) {
        if (verbose) fprintf(stderr, "calc_corr: division by zero \n");
        gamma = 0.0;
    } else {
        gamma = fabs(gamma_num / gamma_denom);
    }

    if (gamma <= 0.5) {
        gamma = 0.0;
    } else {
        gamma = 100.0 * sqrt((gamma * 2.0) - 1.0);
    }

    if (debug)
        fprintf(stdout, " corr %lf \n", gamma);

    return gamma;
}

/*-------------------------------------------------------------------------------*/
void do_time_corr_parallel(struct xcorrp *xc, int iloc, int *i1, int *i2, double *corr) {
    float ipeak = -9999, jpeak = -9999;
    float max_corr = -1;
    const int ysearch = xc->ysearch;
    const int xsearch = xc->xsearch;
    const int nxc = xc->nxc;
    const int nyc = xc->nyc;
    
    /* Parallelize the outer loop with dynamic scheduling */
    #pragma omp parallel for collapse(2) schedule(dynamic) reduction(max:max_corr) \
            private(ipeak, jpeak) if (ysearch > 16 || xsearch > 16)
    for (int ioff = -ysearch; ioff < ysearch; ioff++) {
        for (int joff = -xsearch; joff < xsearch; joff++) {
            const int ic = ioff + ysearch;
            const int jc = joff + xsearch;
            if (ic >= nyc || jc >= nxc) continue;
            double current_corr;

            /* calculate the correlation for each patch */
            if (xc->corr_flag == 0)
                current_corr = calc_time_corr_parallel(xc, ioff, joff, i1, i2);
            else if (xc->corr_flag == 1)
                current_corr = calc_time_corr_hat_parallel(xc, ioff, joff, i1, i2);
            else
                continue;  // Skip invalid corr_flag

            corr[ic * nxc + jc] = current_corr;

            if (current_corr > max_corr) {
                max_corr = current_corr;
                jpeak = joff;
                ipeak = ioff;
            }
        }
    }

    /* Update location data - this is thread-safe as each iloc is unique */
    xc->loc[iloc].xoff = -1 * jpeak;
    xc->loc[iloc].yoff = -1 * ipeak;
    xc->loc[iloc].corr = max_corr;

    if (debug) {
        fprintf(stdout, " (time) jpeak %f xoffset %d \n", jpeak, xc->x_offset);
        fprintf(stdout, " (time) ipeak %f yoffset %d \n", ipeak, xc->y_offset);
    }
}