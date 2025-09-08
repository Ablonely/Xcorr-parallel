# xcorr_parallel

## Description

1. *This project plans to improve the **xcorr function** in GMTSAR using OpenMP.*
2. *However, thread management and memory management in the code are relatively complex, so this function is still under development.*
3. *We hope it can be put into use in the near future.*
## Installation

```shell
# Install fftw
sudo apt-get update
sudo apt-get install libfftw3-dev libfftw3-single3

# Preparation for compilation
# 1. Open the INSTALLATION_PATH/GMTSAR/gmtsar/gmtsar.h
#    Replace content (a) with (b)

# (a)
/* global variables */
int verbose;     /* controls minimal level of output 	*/
int debug;       /* more output 				*/
int swap;        /* whether to swap bytes 		*/
int quad_pol;    /* quad polarization data 		*/
int force_slope; /* whether to force the slope 		*/
int dopp;        /* whether to calculate doppler 	*/
int roi_flag;    /* whether to write roi.in 		*/
int sio_flag;    /* whether to write PRM file 		*/
int nodata;
int quiet_flag;
double forced_slope; /* value to set chirp_slope to		*/
int SAR_mode;        /* 0 => high-res                        */
                     /* 1 => wide obs                        */
                     /* 2 => polarimetry                     */
                     /* from ALOS Product Format 3-2         */
#endif               /* GMTSAR_H */

# (b)
/* global variables */
extern int verbose;
extern int debug;
extern int swap;
extern int quad_pol;
extern int force_slope;
extern int dopp;
extern int roi_flag;
extern int sio_flag;
extern int nodata;
extern int quiet_flag;

extern double forced_slope;
extern int SAR_mode;

# 2. Create a new file (globals.c) and enter the following content:
/* globals.c
 * Define global variables for the GMTSAR project
 * Note: Corresponding to extern declarations in gmtsar.h
 */
#include "gmtsar.h"
/* 全局变量定义 */
int verbose = 0;        /* controls minimal level of output 	*/
int debug = 0;          /* more output 				*/
int swap = 0;           /* whether to swap bytes 		*/
int quad_pol = 0;       /* quad polarization data 		*/
int force_slope = 0;    /* whether to force the slope 		*/
int dopp = 0;           /* whether to calculate doppler 	*/
int roi_flag = 0;       /* whether to write roi.in 		*/
int sio_flag = 0;       /* whether to write PRM file 		*/
int nodata = 0;         
int quiet_flag = 0;    
double forced_slope = 0.0; 
int SAR_mode = 0;  /* SAR mode:
				   /* 0 => high-res                        */
                     /* 1 => wide obs                        */
                     /* 2 => polarimetry                     */
                     /* from ALOS Product Format 3-2         */

# Install xcorr_parallel
git https://github.com/Ablonely/Xcorr-parallel.git
cd /usr/local/GMTSAR/gmtsar/xcorr_para
make

rm /usr/local/GMTSAR/bin/xcorr_parallel
sudo ln -s /usr/local/GMTSAR/gmtsar/xcorr_para/xcorr_parallel /usr/local/GMTSAR/bin/

```
## Usage

```shell
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

# example
xcorr S1Ahh_20170104.PRM S1Ahh_20170110.PRM -time -nx 240 -ny 30 -xsearch 32 -ysearch 32
```


