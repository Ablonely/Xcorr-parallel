/*-------------------------------------------------------*/
/* Parallel version: Calculate pixel locations for each test point */
#include "gmtsar.h"
#include "xcorr_parallel.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void get_locations_parallel(struct xcorrp *xc) {
    int i, j;

    /* Calculate step size (preserve original logic) */
    xc->x_inc = (xc->m_nx - 2 * (xc->xsearch + xc->nx_corr)) / (xc->nxl + 3);
    xc->y_inc = (xc->m_ny - 2 * (xc->ysearch + xc->ny_corr)) / (xc->nyl + 1);

    /* Allocate memory, with one extra column as backup */
    xc->loc = malloc(xc->nyl * (xc->nxl + 1) * sizeof(struct locsp));
    if (!xc->loc) {
        fprintf(stderr, "[PARALLEL] get_locations_parallel: malloc failed\n");
        exit(1);
    }

    /* Total number of locations */
    xc->nlocs = xc->nyl * xc->nxl;

    /* Parallel fill of loc[] array */
    #pragma omp parallel for collapse(2)
    for (j = 1; j <= xc->nyl; j++) {
        for (i = 2; i <= xc->nxl + 1; i++) {
            int idx = (j - 1) * xc->nxl + (i - 2);
            xc->loc[idx].x = xc->npx + i * xc->x_inc;
            xc->loc[idx].y = xc->npy + j * xc->y_inc;
            /* Initialize other fields */
            xc->loc[idx].qflag = 0;
            xc->loc[idx].corr  = 0.0f;
            xc->loc[idx].xoff  = 0.0f;
            xc->loc[idx].yoff  = 0.0f;
            xc->loc[idx].xfrac = 0.0f;
            xc->loc[idx].yfrac = 0.0f;
            xc->loc[idx].m1    = 0;
            xc->loc[idx].m2    = 0;
        }
    }

    fprintf(stderr,
            " locations  n %d nx %d nyl %d nxl %d x_inc %d y_inc %d\n",
            xc->nlocs, xc->m_nx, xc->nyl, xc->nxl,
            xc->x_inc, xc->y_inc);
}