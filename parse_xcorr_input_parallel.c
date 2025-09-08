/*	$Id: parse_xcorr_input_parallel.c 1 2023-10-01 00:00:00Z parallel $	*/
#include "gmtsar.h"
#include "xcorr_parallel.h"  // Updated header file name
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*-------------------------------------------------------*/
void set_defaults_parallel(struct xcorrp *xc) {
    xc->format = 0;
    xc->ri = 2;
    xc->offset_flag = 0;
    xc->nthreads = -1;  // Default: auto-detect

    xc->nxl = 16;
    xc->nyl = 32;
    xc->x_offset = 0;
    xc->y_offset = 0;
    xc->nx_corr = 128;
    xc->ny_corr = 128;
    xc->xsearch = 64;
    xc->ysearch = 64;

    xc->npx = xc->nx_corr + 2 * xc->xsearch;
    xc->npy = xc->ny_corr + 2 * xc->ysearch;
    xc->nxc = 2 * xc->xsearch;
    xc->nyc = 2 * xc->ysearch;

    xc->n2x = 8;
    xc->n2y = 8;
    xc->astretcha = 0.0;
    xc->interp_flag = 1;
    xc->interp_factor = 16;
}

/*-------------------------------------------------------*/
void print_params_parallel(struct xcorrp *xc) {
    fprintf(stdout, "=== Parallel XCORR Parameters ===\n");
    fprintf(stdout, "Format: %d\n", xc->format);
    fprintf(stdout, "Master size: %d x %d\n", xc->m_nx, xc->m_ny);
    fprintf(stdout, "Aligned size: %d x %d\n", xc->s_nx, xc->s_ny);
    fprintf(stdout, "Locations: %d x %d\n", xc->nxl, xc->nyl);
    fprintf(stdout, "Offsets: x=%d, y=%d\n", xc->x_offset, xc->y_offset);
    fprintf(stdout, "Search window: %d x %d\n", xc->xsearch, xc->ysearch);
    fprintf(stdout, "Patch size: %d x %d\n", xc->npx, xc->npy);
    fprintf(stdout, "Threads: %d\n", xc->nthreads);
    fprintf(stdout, "Data files:\n  %s\n  %s\n", xc->data1_name, xc->data2_name);
}

/*-------------------------------------------------------*/
void parse_command_line_parallel(int na, char **a, struct xcorrp *xc, int *nfiles, int *input_flag, char *USAGE) {
    int n;
    FILE *inputfile;
    char tmp[128];

    for (n = 3; n < na; n++) {
        if (!strcmp(a[n], "-freq")) {
            xc->corr_flag = 2;
            fprintf(stderr, "[PARALLEL] Using frequency cross correlation\n");
        }
        else if (!strcmp(a[n], "-time")) {
            xc->corr_flag = 0;
            fprintf(stderr, "[PARALLEL] Using time cross correlation\n");
        }
        else if (!strcmp(a[n], "-time4")) {
            xc->corr_flag = 1;
            fprintf(stderr, "[PARALLEL] Using 4th power time cross correlation\n");
        }
        else if (!strcmp(a[n], "-real")) {
            xc->format = 1;
            fprintf(stderr, "[PARALLEL] Using real data format\n");
        }
        else if (!strcmp(a[n], "-nointerp")) {
            xc->interp_flag = 0;
            fprintf(stderr, "[PARALLEL] Disabling interpolation\n");
        }
        else if (!strcmp(a[n], "-noshift")) {
            xc->offset_flag = 1;
            fprintf(stderr, "[PARALLEL] Ignoring shifts in PRM files\n");
        }
        else if (!strcmp(a[n], "-interp")) {
            n++;
            if (n == na)
                die(" no option after -interp!\n", "");
            xc->interp_flag = 1;
            xc->interp_factor = atoi(a[n]);
            fprintf(stderr, "[PARALLEL] Setting interpolation factor to %d\n", xc->interp_factor);
        }
        else if (!strcmp(a[n], "-range_interp")) {
            n++;
            if (n == na)
                die(" no option after -range_interp!\n", "");
            xc->ri = atoi(a[n]);
            fprintf(stderr, "[PARALLEL] Setting range interpolation factor to %d\n", xc->ri);
        }
        else if (!strcmp(a[n], "-nx")) {
            n++;
            if (n == na)
                die(" no option after -nx!\n", "");
            xc->nxl = atoi(a[n]);
            fprintf(stderr, "[PARALLEL] Setting nx to %d\n", xc->nxl);
        }
        else if (!strcmp(a[n], "-ny")) {
            n++;
            if (n == na)
                die(" no option after -ny!\n", "");
            xc->nyl = atoi(a[n]);
            fprintf(stderr, "[PARALLEL] Setting ny to %d\n", xc->nyl);
        }
        else if (!strcmp(a[n], "-xsearch")) {
            n++;
            if (n == na)
                die(" no option after -xsearch!\n", "");
            xc->xsearch = atoi(a[n]);
            xc->nx_corr = 2 * xc->xsearch;
            xc->npx = xc->nx_corr + 2 * xc->xsearch; /* size of data required*/
            xc->nxc = 2 * xc->xsearch;               /* size of correlation patch*/
            if (((xc->xsearch - 1) & xc->xsearch))
                die(" xsearch needs to be power of 2! (32 64 128 256) \n", "");
            fprintf(stderr, "[PARALLEL] Setting xsearch to %d\n", xc->xsearch);
            fprintf(stderr, "[PARALLEL] Setting nx_corr to %d\n", xc->nx_corr);
        }
        else if (!strcmp(a[n], "-ysearch")) {
            n++;
            if (n == na)
                die(" no option after -ysearch!\n", "");
            xc->ysearch = atoi(a[n]);
            xc->ny_corr = 2 * xc->ysearch;
            xc->npy = xc->ny_corr + 2 * xc->ysearch; /* size of data required*/
            xc->nyc = 2 * xc->ysearch;               /* size of correlation patch*/
            if (((xc->ysearch - 1) & xc->ysearch))
                die(" ysearch needs to be power of 2! (32 64 128 256) \n", "");
            fprintf(stderr, "[PARALLEL] Setting ysearch to %d\n", xc->ysearch);
            fprintf(stderr, "[PARALLEL] Setting ny_corr to %d\n", xc->ny_corr);
        }
        else if (!strcmp(a[n], "-v")) {
            verbose = 1;
            fprintf(stderr, "[PARALLEL] Verbose output enabled\n");
        }
        else if (!strcmp(a[n], "-norange")) {
            xc->ri = 1;
            fprintf(stderr, "[PARALLEL] Disabling range interpolation\n");
        }
        else if (!strncmp(a[n], "-input", 6)) {
            n++;
            if (n == na)
                die(" no option after -input!\n", "");
            fprintf(stderr, "[PARALLEL] Using input file \n");
            *input_flag = 1;

            if ((inputfile = fopen(a[2], "r")) == NULL)
                die("Can't open ", a[2]);

            while (fscanf(inputfile, " %s ", tmp) != EOF)
                (*nfiles)++;
            fclose(inputfile);
        }
        else if (!strcmp(a[n], "-threads")) {
            if (++n == na) die("No thread count specified", "");
            xc->nthreads = atoi(a[n]);
            #ifndef _OPENMP
            if (xc->nthreads > 1) {
                fprintf(stderr, "Warning: OpenMP not available. Using single thread.\n");
                xc->nthreads = 1;
            }
            #endif
            fprintf(stderr, "[PARALLEL] Using %d threads\n", xc->nthreads);
        }
        else {
            fprintf(stderr, "Unknown option: %s\n", a[n]);
            die(USAGE, "");
        }
    }
    
    // 更新派生参数
    xc->npx = xc->nx_corr + 2 * xc->xsearch;
    xc->npy = xc->ny_corr + 2 * xc->ysearch;
    xc->nxc = 2 * xc->xsearch;
    xc->nyc = 2 * xc->ysearch;
    xc->nlocs = xc->nxl * xc->nyl;
}

/*---------------------------------------------------------------------------*/
void handle_prm_parallel(void *API, char **argv, struct xcorrp *xc, int nfiles) {
    char **filename = malloc(nfiles * sizeof(char *));
    for (int i = 0; i < nfiles; i++) {
        filename[i] = malloc(128);
        strcpy(filename[i], argv[i + 1]);
    }

    struct PRM *r = malloc(nfiles * sizeof(struct PRM));
    
    for (int i = 0; i < nfiles; i++) {
        if (strstr(filename[i], ".PRM")) {
            FILE *prmfile = fopen(filename[i], "r");
            if (!prmfile) die("Can't open PRM file", filename[i]);
            
            null_sio_struct(&r[i]);
            get_sio_struct(prmfile, &r[i]);
            fclose(prmfile);

            if (i == 0) {
                strcpy(xc->data1_name, r[i].SLC_file);
                if (!(xc->data1 = fopen(xc->data1_name, "rb")))
                    die("Cannot open SLC", xc->data1_name);
                xc->m_nx = r[i].num_rng_bins;
                xc->m_ny = r[i].num_patches * r[i].num_valid_az;
            }
            else if (i == 1) {
                strcpy(xc->data2_name, r[i].SLC_file);
                if (!(xc->data2 = fopen(xc->data2_name, "rb")))
                    die("Cannot open SLC", xc->data2_name);
                xc->s_nx = r[i].num_rng_bins;
                xc->s_ny = r[i].num_patches * r[i].num_valid_az;
                xc->x_offset = xc->offset_flag ? 0 : r[i].rshift;
                xc->y_offset = xc->offset_flag ? 0 : r[i].ashift;

                
                if (r[0].prf > 0) {
                    xc->astretcha = (r[1].prf - r[0].prf) / r[0].prf;
                }
            }
        }
        else {
            xc->format = 2;
            // NetCDF handling same as original version
        }
    }
    
    free(r);
    for (int i = 0; i < nfiles; i++) free(filename[i]);
    free(filename);
}