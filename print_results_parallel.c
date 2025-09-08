#include "gmtsar.h"
#include "siocomplex.h"
#include "xcorr_parallel.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Thread-safe version: Print results for one location */
void print_results_parallel(struct xcorrp *xc, int iloc) {
    int ishft;
    int interp;
    float xoff, xfrac;
    float yoff, yfrac;
    float corr;

    xoff  = xc->loc[iloc].xoff;
    xfrac = xc->loc[iloc].xfrac;
    yoff  = xc->loc[iloc].yoff;
    yfrac = xc->loc[iloc].yfrac;
    corr  = xc->loc[iloc].corr;

    if (xc->interp_flag)
        interp = xc->interp_factor;
    else
        interp = 1;

    ishft = (int)xc->loc[iloc].y * xc->astretcha;

    if (debug)
        fprintf(stdout,
                " xoff %f xfrac %f xoff/ri %f rshift %d "
                "yoff %f yfrac %f ashift %d\n",
                xoff, xfrac, xoff / (float)xc->ri,
                xc->x_offset, yoff, yfrac, xc->y_offset);

    /* Apply offset corrections */
    xoff = (xoff / (float)xc->ri) - (xfrac / (float)xc->ri) + xc->x_offset;
    yoff = yoff - yfrac + xc->y_offset + ishft;

    if (verbose) {
        fprintf(stdout,
                " location %d (%3d,%3d) interpolation (range %d corr %d) "
                "correlation %6.2f offset (%6.3f,%6.3f)\n",
                iloc, xc->loc[iloc].x, xc->loc[iloc].y,
                xc->ri, interp, corr, xoff, yoff);
    }

    /* File writing section with lock to avoid parallel conflicts */
    #pragma omp critical(print_results_file)
    {
        fprintf(xc->file, " %d %6.3f %d %6.3f %6.2f \n",
                xc->loc[iloc].x, xoff,
                xc->loc[iloc].y, yoff, corr);
    }
}

/* Print complex matrix to stdout (thread-safe but I/O serialized) */
void print_complex_parallel(struct FCOMPLEX *a, int ny, int nx, int real_flag) {
    #pragma omp critical(print_complex)
    {
        if (real_flag == 0)
            fprintf(stdout, "\ncomplex: \n");
        if (real_flag == 1)
            fprintf(stdout, "\ncomplex (real only): \n");

        for (int i = 0; i < ny; i++) {
            for (int j = 0; j < nx; j++) {
                if (real_flag == 0)
                    fprintf(stdout, "(%6.2f,%6.2f) ",
                            a[i * nx + j].r, a[i * nx + j].i);
                else
                    fprintf(stdout, "%3.1f ", a[i * nx + j].r);
            }
            fprintf(stdout, "\n");
        }
        fprintf(stdout, "\n");
    }
}

/*-------------------------------------------------------------------------------*/
void print_float_parallel(float *a, int ny, int nx) {
    #pragma omp critical(print_float)
    {
        for (int i = 0; i < ny; i++) {
            for (int j = 0; j < nx; j++)
                fprintf(stdout, " %4.2f ", a[i * nx + j]);
            fprintf(stdout, "\n");
        }
    }
}

/*-------------------------------------------------------------------------------*/
void print_double_parallel(double *a, int ny, int nx) {
    #pragma omp critical(print_double)
    {
        for (int i = 0; i < ny; i++) {
            for (int j = 0; j < nx; j++)
                fprintf(stdout, " %4.2lf ", a[i * nx + j]);
            fprintf(stdout, "\n");
        }
    }
}

/*-------------------------------------------------------------------------------*/
void print_int_parallel(int *a, int ny, int nx) {
    #pragma omp critical(print_int)
    {
        for (int i = 0; i < ny; i++) {
            for (int j = 0; j < nx; j++)
                fprintf(stdout, " %d ", a[i * nx + j]);
            fprintf(stdout, "\n");
        }
    }
}