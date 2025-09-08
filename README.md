# xcorr_parallel

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
git xxx
cd /usr/local/GMTSAR/gmtsar/xcorr_para
make

rm /usr/local/GMTSAR/bin/xcorr_parallel
sudo ln -s /usr/local/GMTSAR/gmtsar/xcorr_para/xcorr_parallel /usr/local/GMTSAR/bin/

```

