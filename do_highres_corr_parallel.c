/*-------------------------------------------------------*/
/* Parallel version: Sub-pixel high-resolution correlation (FFT interpolation method) */
#include "gmtsar.h"
#include "lib_functions.h"
#include "xcorr_parallel.h"
#include <math.h>
#include <string.h>
#include <omp.h>

void do_highres_corr_parallel(void *API, struct xcorrp *xc, int iloc, void *c3_unused /* Keep consistent with header */) {
    (void)c3_unused;  /* Unused parameter, prevent warning */

    int i, j, ic, jc, k;
    int nx, ny, nx2, ny2, ifc;

    float max_corr = -1.0f;
    float ipeak = -9999.0f, jpeak = -9999.0f;
    float sub_xoff, sub_yoff;

    ifc = xc->interp_factor;

    /* Complex correlation window dimensions (must be power of 2) */
    nx = xc->n2x;
    ny = xc->n2y;

    /* Matrix dimensions after interpolation */
    nx2 = ifc * nx;
    ny2 = ifc * ny;

    /* Center window around 1-pixel resolution offset (ri factor causes x-direction to have factor of 2) */
    jc = (xc->nxc / 2) - nx / 2 - (int)xc->loc[iloc].xoff;
    ic = (xc->nyc / 2) - ny / 2 - (int)xc->loc[iloc].yoff;

    /* Thread-local buffers to avoid shared data race conditions */
    struct FCOMPLEX *md_local      = (struct FCOMPLEX *)malloc((size_t)nx * ny   * sizeof(struct FCOMPLEX));
    struct FCOMPLEX *cd_exp_local  = (struct FCOMPLEX *)malloc((size_t)nx2 * ny2 * sizeof(struct FCOMPLEX));
    if (!md_local || !cd_exp_local) {
        fprintf(stderr, "[PARALLEL] do_highres_corr_parallel: malloc failed (nx=%d ny=%d nx2=%d ny2=%d)\n",
                nx, ny, nx2, ny2);
        free(md_local);
        free(cd_exp_local);
        return;
    }

    /* Extract window from correlation matrix (double) to complex buffer, apply 4th root to amplitude to enhance peak resolution */
    /* Note: Boundary check, skip if out of bounds (keep as 0) */
    #pragma omp parallel for collapse(2) if (nx*ny >= 4096)
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            struct FCOMPLEX val = {0.0f, 0.0f};
            int rr = ic + i;
            int cc = jc + j;
            if (rr >= 0 && rr < xc->nyc && cc >= 0 && cc < xc->nxc) {
                k = rr * xc->nxc + cc;
                /* Original implementation: powf(corr, 0.25) in real part, imaginary part as 0 */
                float amp = (float)pow(xc->corr[k], 0.25);
                val.r = amp;
                val.i = 0.0f;
            }
            md_local[i * nx + j] = val;
        }
    }

    if (debug) {
        print_complex(md_local, nx, ny, 1);
    }

    /* 2D frequency domain interpolation: interpolate md_local from (ny,nx) to (ny2,nx2) -> cd_exp_local */
    #ifdef XCORR_INTERP_NEEDS_LOCK
    #pragma omp critical(gmt_fft_interp)
    {
        fft_interpolate_2d(API, md_local, ny, nx, cd_exp_local, ny2, nx2, ifc);
    }
    #else
    fft_interpolate_2d(API, md_local, ny, nx, cd_exp_local, ny2, nx2, ifc);
    #endif

    if (ifc <= 4 && debug) {
        print_complex(cd_exp_local, nx2, ny2, 1);
    }

    /* Search for maximum value and its position in the interpolated matrix (relative to center-aligned coordinate system) */
    /* Use per-thread local maximum aggregation to avoid frequent critical section entry */
    #pragma omp parallel
    {
        float t_max = -1.0f;
        int   t_i   = -1;
        int   t_j   = -1;

        #pragma omp for nowait collapse(2)
        for (i = 0; i < ny2; i++) {
            for (j = 0; j < nx2; j++) {
                float r = cd_exp_local[i * nx2 + j].r;
                if (r > t_max) {
                    t_max = r;
                    t_i = i;
                    t_j = j;
                }
            }
        }

        #pragma omp critical(update_peak)
        {
            if (t_max > max_corr) {
                max_corr = t_max;
                /* Note: Original code shifts i/j relative to center here */
                ipeak = (float)(t_i - (ny2 / 2));
                jpeak = (float)(t_j - (nx2 / 2)); 
            }
        }
    }

    /* Calculate sub-pixel offset (in pixel units) */
    sub_xoff = jpeak / (float)ifc;
    sub_yoff = ipeak / (float)ifc;

    if (debug) {
        fprintf(stderr,
                " highres [ri %d ifc %d](%4.1f %4.1f) (nx %d ny %d -> nx2 %d ny2 %d) "
                "jpeak %f ipeak %f : %f %f : %4.2f\n",
                xc->ri, ifc, xc->loc[iloc].xoff, xc->loc[iloc].yoff,
                nx, ny, nx2, ny2, jpeak, ipeak, sub_xoff, sub_yoff, max_corr);
    }

    /* Write back sub-pixel results */
    xc->loc[iloc].xfrac = sub_xoff;
    xc->loc[iloc].yfrac = sub_yoff;

    /* Free thread-local buffers */
    free(md_local);
    free(cd_exp_local);
}