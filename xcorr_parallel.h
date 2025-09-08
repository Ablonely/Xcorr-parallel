/*	$Id: xcorr_parallel.h 39 2013-04-07 00:49:34Z pwessel $	*/
#ifndef XCORR_PARALLEL_H
#define XCORR_PARALLEL_H

#include <stdio.h>
#include <omp.h>  // Add OpenMP support

// Thread-local data structure
struct thread_data {
    short *mask;             // Per-thread independent mask array
    int *i1;                 // Per-thread independent data matrix 1
    int *i2;                 // Per-thread independent data matrix 2
    struct FCOMPLEX *c1;     // Per-thread independent data block 1
    struct FCOMPLEX *c2;     // Per-thread independent data block 2
    struct FCOMPLEX *c3;     // Per-thread independent cross-correlation result
};

struct locsp {
    int x;       /* x pixel location */
    int y;       /* y pixel location */
    int qflag;   /* quality flag 1 = good */
    float corr;  /* correlation, normalized time domain */
    float yoff;  /* estimated y offset */
    float xoff;  /* estimated x offset */
    float xfrac; /* estimated x offset (fraction) */
    float yfrac; /* estimated y offset (fraction) */
    int m1;      /* mean value */
    int m2;      /* mean value */
};

struct xcorrp {
    // Parallel computation related members
    int nthreads;            /* Number of threads: >0 user specified, -1 auto-detected */
    struct thread_data *thread_data; /* Thread-local storage array */
    
    // Original members remain unchanged
    int format;              /* type of input data [0 short complex, 1 real float, 3 real grd] */
    int corr_flag;           /* type of correlation flag 0 = standard; 1 = Gatelli; 2 = fft */
    int offset_flag;         /* set offset to zero (ignore prm)  */
    int interp_flag;         /* interpolation flag 1 = yes 0 = no */
    int interp_factor;       /* interpolation factor (power of 2) */
    int nx_corr;             /* size of correlation window */
    int ny_corr;             /* size of correlation window */
    int xsearch;             /* size of x search offset */
    int ysearch;             /* size of y search offset */
    int m_nx;                /* x size of master file */
    int m_ny;                /* y size of master file */
    int s_nx;                /* x size of aligned file */
    int s_ny;                /* y size of aligned file */
    int nxc;                 /* x size of correlation */
    int nyc;                 /* y size of correlation */
    int npx;                 /* x size of patch */
    int npy;                 /* y size of patch */
    int nxl;                 /* x number of locations */
    int nyl;                 /* y number of locations */
    int x_offset;            /* intial starting offset in x */
    int y_offset;            /* intial starting offset in y */
    int nlocs;               /* number of locations */
    int x_inc;               /* x distance between locations */
    int y_inc;               /* y distance between locations */
    int ri;                  /* range interpolation factor (must be power of two) */
    short *mask;             /* master mask file (short integer) */
    int *i1;                 /* Main data matrix 1 (integer) */
    int *i2;                 /* Main data matrix 2 (integer) */
    int n2x;                 /* size of interpolation */
    int n2y;                 /* size of interpolation */
    struct FCOMPLEX *d1;     /* Main data 1 (amplitude in real, imag = 0)*/
    struct FCOMPLEX *d2;     /* Main data 2 (amplitude in real, imag = 0)*/
    struct FCOMPLEX *c1;     /* Main data block 1 (complex float) */
    struct FCOMPLEX *c2;     /* Main data block 2 (complex float) */
    struct FCOMPLEX *c3;     /* Main cross-correlation result (complex float) */
    struct FCOMPLEX *ritmp;  /* tmp array for range interpoation */
    double *corr;            /* correlation (real double) */
    double astretcha;        /* azimuth stretch parameter estimated from (prf2-prf1)/prf1 */
    struct FCOMPLEX *md;     /* interpolation file */
    struct FCOMPLEX *md_exp; /* interpolation file */
    struct FCOMPLEX *sd;     /* interpolation file */
    struct FCOMPLEX *sd_exp; /* interpolation file */
    struct FCOMPLEX *cd_exp; /* interpolation file */
    double *interp_corr;     /* interpolation file */
    FILE *data1;             /* data file 1 */
    FILE *data2;             /* data file 2 */
    FILE *param;             /* input parameters file */
    FILE *file;              /* output file (offsets) */
    char param_name[128];
    char data1_name[128];
    char data2_name[128];
    char filename[128];
    struct locsp *loc;
    struct GMT_GRID *D1;
    struct GMT_GRID *D2;
};


void allocate_arrays_parallel(struct xcorrp *xc);
void make_mask_parallel(struct xcorrp *xc, short *mask);
void do_time_corr_parallel(struct xcorrp *xc, int iloc, int *i1, int *i2, double *corr);
void do_freq_corr_parallel(void *API, struct xcorrp *xc, int iloc);
void do_highres_corr_parallel(void *API, struct xcorrp *xc, int iloc, void *c3_unused);
void print_results_parallel(struct xcorrp *xc, int iloc);
void set_defaults_parallel(struct xcorrp *xc);
void parse_command_line_parallel(int argc, char **argv, struct xcorrp *xc, int *nsta, int *nch, char *usage);
void handle_prm_parallel(void *API, char **argv, struct xcorrp *xc, int format);
void print_params_parallel(struct xcorrp *xc);
void get_locations_parallel(struct xcorrp *xc);
void read_xcorr_data_parallel(struct xcorrp *xc, int iloc);
void assign_values_parallel(void *API, struct xcorrp *xc, int iloc, 
                   struct FCOMPLEX *c1, struct FCOMPLEX *c2, 
                   struct FCOMPLEX *ritmp, int *i1, int *i2, short *mask);
int do_range_interpolate_parallel(void *API, struct FCOMPLEX *c, int nx, int ri, struct FCOMPLEX *work);
double calc_time_corr_parallel(struct xcorrp *xc, int ioff, int joff, int *i1, int *i2);

// Declare functions in the GMTSAR library
void fft_interpolate_1d(void *API, struct FCOMPLEX *c, int nx, struct FCOMPLEX *work, int ri);
float Cabs(struct FCOMPLEX z);
struct FCOMPLEX Cmul(struct FCOMPLEX z1, struct FCOMPLEX z2);
struct FCOMPLEX Conjg(struct FCOMPLEX z);
struct FCOMPLEX RCmul(float d, struct FCOMPLEX z);
void die(char *message, char *arg);
void print_complex(struct FCOMPLEX *c, int nx, int ny, int real_flag);

#endif /* XCORR_PARALLEL_H */