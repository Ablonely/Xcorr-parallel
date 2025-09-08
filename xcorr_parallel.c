/*	$Id: xcorr_parallel.c 1 2023-10-01 00:00:00Z parallel $	*/
/***************************************************************************/
/* Parallel version of xcorr with OpenMP optimization                      */
/***************************************************************************/

#include "gmtsar.h"
#include "xcorr_parallel.h"  // Updated header file name
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif


// Add missing function declarations
//void do_time_corr_parallel(struct xcorr *xc, int iloc, int *i1, int *i2, struct FCOMPLEX *c3);
//void do_freq_corr_parallel(void *API, struct xcorr *xc, int iloc, struct FCOMPLEX *c1, struct FCOMPLEX *c2, struct FCOMPLEX *c3);
//void do_highres_corr_parallel(void *API, struct xcorr *xc, int iloc, struct FCOMPLEX *c3);
//void print_results_parallel(struct xcorr *xc, int iloc);
//void get_locations_parallel(struct xcorr *xc);

char *USAGE = "xcorr_parallel [GMTSAR] - Parallel Compute 2-D cross-correlation of two images\n\n"
              "\nUsage: xcorr master.PRM aligned.PRM [-time] [-real] [-freq] [-nx n] [-ny n]"
              "[-xsearch xs] [-ysearch ys] [-threads n]\n"
              "master.PRM     	PRM file for reference image\n"
              "aligned.PRM      PRM file of secondary image\n"
              "-time      		use time cross-correlation\n"
              "-freq      		use frequency cross-correlation (default)\n"
              "-real      		read float numbers instead of complex numbers\n"
              "-noshift  		ignore ashift and rshift in prm file (set to 0)\n"
              "-nx  nx    		number of locations in x (range) direction (int)\n"
              "-ny  ny    		number of locations in y (azimuth) direction (int)\n"
              "-nointerp     		do not interpolate correlation function\n"
              "-range_interp ri  	interpolate range by ri (power of two) [default: 2]\n"
              "-norange     		do not range interpolate \n"
              "-xsearch xs		search window size in x (range) direction (int "
              "power of 2 [32 64 128 256])\n"
              "-ysearch ys		search window size in y (azimuth) direction "
              "(int power of 2 [32 64 128 256])\n"
              "-interp  factor    	interpolate correlation function by factor "
              "(int) [default, 16]\n"
              "-v			verbose\n"
              "output: \n freq_xcorr.dat (default) \n time_xcorr.dat (if -time option))\n"
              "\nuse fitoffset.csh to convert output to PRM format\n"
              "\nExample:\n"
              "xcorr IMG-HH-ALPSRP075880660-H1.0__A.PRM "
              "IMG-HH-ALPSRP129560660-H1.0__A.PRM -nx 20 -ny 50 \n"
              "xcorr file1.grd file2.grd -nx 20 -ny 50 (takes grids with real numbers)\n";

/*-------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/
int do_range_interpolate_parallel(void *API, struct FCOMPLEX *c, int nx, int ri, struct FCOMPLEX *work) {
    int i;
    fft_interpolate_1d(API, c, nx, work, ri);
    for (i = 0; i < nx; i++) {
        c[i].r = work[i + nx / 2].r;
        c[i].i = work[i + nx / 2].i;
    }
    return (EXIT_SUCCESS);
}

/*-------------------------------------------------------------------------------*/
void assign_values_parallel(void *API, struct xcorrp *xc, int iloc, 
                   struct FCOMPLEX *c1, struct FCOMPLEX *c2, 
                   struct FCOMPLEX *ritmp, int *i1, int *i2, short *mask) {
    int i, j, k, sx, mx;
    double mean1 = 0.0, mean2 = 0.0;

    mx = xc->loc[iloc].x - xc->npx / 2;
    sx = xc->loc[iloc].x + xc->x_offset - xc->npx / 2;

    // Optimize loop structure
    #pragma omp simd collapse(2)
    for (i = 0; i < xc->npy; i++) {
        for (j = 0; j < xc->npx; j++) {
            k = i * xc->npx + j;
            c1[k] = xc->d1[i * xc->m_nx + mx + j];
            c2[k] = xc->d2[i * xc->s_nx + sx + j];
        }
    }

    if (xc->ri > 1) {
        #pragma omp parallel for
        for (i = 0; i < xc->npy; i++) {
            do_range_interpolate_parallel(API, &c1[i * xc->npx], xc->npx, xc->ri, ritmp);
            do_range_interpolate_parallel(API, &c2[i * xc->npx], xc->npx, xc->ri, ritmp);
        }
    }

    // Vectorized processing
    #pragma omp simd reduction(+:mean1, mean2)
    for (i = 0; i < xc->npy * xc->npx; i++) {
        c1[i].r = Cabs(c1[i]);
        c1[i].i = 0.0f;
        c2[i].r = Cabs(c2[i]);
        c2[i].i = 0.0f;
        mean1 += c1[i].r;
        mean2 += c2[i].r;
    }

    mean1 /= (xc->npy * xc->npx);
    mean2 /= (xc->npy * xc->npx);

    #pragma omp simd
    for (i = 0; i < xc->npy * xc->npx; i++) {
        c1[i].r -= (float)mean1;
        c2[i].r -= (float)mean2;
        c1[i].i = c2[i].i = 0.0f;
        c2[i].r *= (float)mask[i];
        i1[i] = (int)(c1[i].r);
        i2[i] = (int)(c2[i].r);
    }
}

/*-------------------------------------------------------------------------------*/
void read_xcorr_data_parallel(struct xcorrp *xc, int iloc) {
    size_t block_size = xc->npy * xc->m_nx * sizeof(struct FCOMPLEX);
    off_t offset_m = (off_t)(xc->loc[iloc].y - xc->npy/2) * xc->m_nx * sizeof(struct FCOMPLEX);
    
    /* Block read master image */
    fseek(xc->data1, offset_m, SEEK_SET);
    if (fread(xc->d1, 1, block_size, xc->data1) != block_size && !feof(xc->data1)) {
        die("Read error (master)", xc->data1_name);
    }
    
    /* Calculate secondary image offset */
    int ishft = (int)(xc->loc[iloc].y * xc->astretcha);
    int y_aligned = xc->loc[iloc].y + xc->y_offset + ishft - xc->npy/2;
    off_t offset_s = (off_t)y_aligned * xc->s_nx * sizeof(struct FCOMPLEX);
    
    /* Block read secondary image */
    fseek(xc->data2, offset_s, SEEK_SET);
    if (fread(xc->d2, 1, block_size, xc->data2) != block_size && !feof(xc->data2)) {
        die("Read error (aligned)", xc->data2_name);
    }
    
    /* Format conversion */
    if (xc->format == 0) {
        #pragma omp parallel for
        for (int i = 0; i < xc->npy * xc->m_nx; i++) {
            xc->d1[i].r = (float)((short *)xc->d1)[2*i];
            xc->d1[i].i = (float)((short *)xc->d1)[2*i+1];
        }
        #pragma omp parallel for
        for (int i = 0; i < xc->npy * xc->s_nx; i++) {
            xc->d2[i].r = (float)((short *)xc->d2)[2*i];
            xc->d2[i].i = (float)((short *)xc->d2)[2*i+1];
        }
    }

    else if (xc->format == 1) {
        #pragma omp parallel for
        for (int i = 0; i < xc->npy * xc->m_nx; i++) {
            xc->d1[i].r = ((float *)xc->d1)[i];
            xc->d1[i].i = 0.0f;
        }
        #pragma omp parallel for
        for (int i = 0; i < xc->npy * xc->s_nx; i++) {
            xc->d2[i].r = ((float *)xc->d2)[i];
            xc->d2[i].i = 0.0f;
        }
    }
}

/*-------------------------------------------------------------------------------*/
void do_correlation_parallel(void *API, struct xcorrp *xc) {
    int total_locations = xc->nyl * xc->nxl;
    allocate_arrays_parallel(xc);

    #ifdef _OPENMP
    int nthreads = xc->nthreads > 0 ? xc->nthreads : (omp_get_num_procs() - 2);
    if (nthreads < 1) nthreads = 1;
    omp_set_num_threads(nthreads);
    #else
    int nthreads = 1;
    #endif
    
    if (verbose) fprintf(stderr, "Processing %d locations with %d threads\n", total_locations, nthreads);
    
    struct ThreadData {
        struct FCOMPLEX *c1, *c2, *c3, *ritmp;
        short *mask;
        int *i1, *i2;
    } *thread_data = malloc(nthreads * sizeof(struct ThreadData));
    
    #pragma omp parallel num_threads(nthreads)
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        
        thread_data[tid].c1 = malloc(xc->npx * xc->npy * sizeof(struct FCOMPLEX));
        thread_data[tid].c2 = malloc(xc->npx * xc->npy * sizeof(struct FCOMPLEX));
        thread_data[tid].c3 = malloc(xc->npx * xc->npy * sizeof(struct FCOMPLEX));
        thread_data[tid].ritmp = malloc(xc->ri * xc->npx * sizeof(struct FCOMPLEX));
        thread_data[tid].mask = malloc(xc->npx * xc->npy * sizeof(short));
        thread_data[tid].i1 = malloc(xc->npx * xc->npy * sizeof(int));
        thread_data[tid].i2 = malloc(xc->npx * xc->npy * sizeof(int));
        
        make_mask_parallel(xc, thread_data[tid].mask);
    }
    
    double start_time = omp_get_wtime();
    
    #pragma omp parallel for num_threads(nthreads) schedule(dynamic)

    for (int iloc = 0; iloc < total_locations; iloc++) {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        
        struct ThreadData *td = &thread_data[tid];
        
        #pragma omp critical
        read_xcorr_data_parallel(xc, iloc);
        
        assign_values_parallel(API, xc, iloc, td->c1, td->c2, td->ritmp, td->i1, td->i2, td->mask);
        
        if (xc->corr_flag < 2)
            //do_time_corr_parallel(xc, iloc, td->i1, td->i2, td->c3);
            do_time_corr_parallel(xc, iloc, td->i1, td->i2, xc->corr);
        if (xc->corr_flag == 2)
            //do_freq_corr_parallel(API, xc, iloc, td->c1, td->c2, td->c3);
            do_freq_corr_parallel(API, xc, iloc);
        if (xc->interp_flag == 1)
            do_highres_corr_parallel(API, xc, iloc, td->c3);
        
        #pragma omp critical
        print_results_parallel(xc, iloc);
    }
    
    double end_time = omp_get_wtime();
    if (verbose) fprintf(stderr, "Parallel time: %.2f seconds\n", end_time - start_time);
    
    #pragma omp parallel num_threads(nthreads)
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        free(thread_data[tid].c1);
        free(thread_data[tid].c2);
        free(thread_data[tid].c3);
        free(thread_data[tid].ritmp);
        free(thread_data[tid].mask);
        free(thread_data[tid].i1);
        free(thread_data[tid].i2);
    }
    free(thread_data);
}

/*-------------------------------------------------------------------------------*/
void make_mask_parallel(struct xcorrp *xc, short *mask) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < xc->npy; i++) {
        for (int j = 0; j < xc->npx; j++) {
            int idx = i * xc->npx + j;
            mask[idx] = 1;
            if (i < xc->ysearch || i >= (xc->npy - xc->ysearch) || 
                j < xc->xsearch || j >= (xc->npx - xc->xsearch)) {
                mask[idx] = 0;
            }
        }
    }
}

/*-------------------------------------------------------------------------------*/
void allocate_arrays_parallel(struct xcorrp *xc) {
    #ifdef _OPENMP
    if (posix_memalign((void**)&xc->d1, 64, xc->m_nx * xc->npy * sizeof(struct FCOMPLEX)) != 0) {
        die("Memory alignment failed for d1", "");
    }
    if (posix_memalign((void**)&xc->d2, 64, xc->s_nx * xc->npy * sizeof(struct FCOMPLEX)) != 0) {
        die("Memory alignment failed for d2", "");
    }
    #else
    xc->d1 = malloc(xc->m_nx * xc->npy * sizeof(struct FCOMPLEX));
    xc->d2 = malloc(xc->s_nx * xc->npy * sizeof(struct FCOMPLEX));
    #endif

    xc->c1 = malloc(xc->npx * xc->npy * sizeof(struct FCOMPLEX));
    xc->c2 = malloc(xc->npx * xc->npy * sizeof(struct FCOMPLEX));
    xc->c3 = malloc(xc->npx * xc->npy * sizeof(struct FCOMPLEX));
    xc->ritmp = malloc(xc->ri * xc->npx * sizeof(struct FCOMPLEX));
    xc->mask = malloc(xc->npx * xc->npy * sizeof(short));
    xc->corr = malloc(2 * xc->ri * xc->nxc * xc->nyc * sizeof(double));

    if (xc->interp_flag == 1) {
        int nx = 2 * xc->n2x;
        int ny = 2 * xc->n2y;
        xc->md = malloc(nx * ny * sizeof(struct FCOMPLEX));
        xc->cd_exp = malloc(nx * xc->interp_factor * ny * xc->interp_factor * sizeof(struct FCOMPLEX));
    }
}

/*-------------------------------------------------------*/
int main(int argc, char **argv) {
    struct xcorrp *xc = calloc(1, sizeof(struct xcorrp));
    void *API = NULL;
    
    /* Timing variables */
    clock_t start, end;
    double cpu_time;
    #ifdef _OPENMP
    double wall_start, wall_end, wall_time;
    #endif

    if ((API = GMT_Create_Session(argv[0], 0U, 0U, NULL)) == NULL)
        return EXIT_FAILURE;

    if (argc < 3) die(USAGE, "");

    verbose = 0;
    debug = 0;
    xc->interp_flag = 0;
    xc->corr_flag = 2;
    xc->nthreads = -1;

    set_defaults_parallel(xc);
    parse_command_line_parallel(argc, argv, xc, &(int){2}, &(int){0}, USAGE);

    if (xc->format != 2) handle_prm_parallel(API, argv, xc, 2);

    if (debug) print_params_parallel(xc);

    strcpy(xc->filename, (xc->corr_flag == 2) ? "freq_xcorr.dat" : 
                         (xc->corr_flag == 1) ? "time_xcorr_Gatelli.dat" : "time_xcorr.dat");
    
    if (!(xc->file = fopen(xc->filename, "w")))
        die("Can't open output file", xc->filename);

    fprintf(stderr, "----------------------------------\n");
    fprintf(stderr, "Result file created successfully!\n");
    
    get_locations_parallel(xc);

    fprintf(stderr, "----------------------------------\n");
    fprintf(stderr, "Monitoring point location successfully identified!\n");

    /* Start timing */
    start = clock();
    #ifdef _OPENMP
    wall_start = omp_get_wtime();
    #endif

    do_correlation_parallel(API, xc);

    fprintf(stderr, "----------------------------------\n");
    fprintf(stderr, "Cross-correlation calculation completed!\n");

    /* End timing and print results */
    end = clock();
    cpu_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    #ifdef _OPENMP
    wall_end = omp_get_wtime();
    wall_time = wall_end - wall_start;
    fprintf(stdout, "CPU time: %.2f seconds\n", cpu_time);
    fprintf(stdout, "Wall-clock time: %.2f seconds\n", wall_time);
    fprintf(stdout, "Parallel efficiency: %.1f%%\n", 
            (cpu_time / (xc->nthreads > 0 ? xc->nthreads : omp_get_num_procs())) / wall_time * 100);
    #else
    fprintf(stdout, "Elapsed time: %.2f seconds\n", cpu_time);
    #endif

    if (xc->format == 0 || xc->format == 1) {
        fclose(xc->data1);
        fclose(xc->data2);
    }

    if (GMT_Destroy_Session(API)) return EXIT_FAILURE;
    return EXIT_SUCCESS;
}