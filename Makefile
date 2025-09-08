# Makefile for xcorr_parallel in xcorr_para/
CC      = gcc
CFLAGS  = -O2 -Wall -fopenmp -I../ -I/usr/include -I/usr/include/gmt
LDFLAGS = -L../ -L/usr/lib/x86_64-linux-gnu
LIBS    = -lgmtsar -lgmt -lm -ltiff -lfftw3f -fopenmp

# Object files
OBJS = xcorr_parallel.o \
       do_freq_xcorr_parallel.o \
       do_highres_corr_parallel.o \
       do_time_int_xcorr_parallel.o \
       get_locations_parallel.o \
       parse_xcorr_input_parallel.o \
       print_results_parallel.o \
       ../globals.o

# Target
TARGET = xcorr_parallel

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

# Compilation rules
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

../globals.o: ../globals.c
	$(CC) $(CFLAGS) -c ../globals.c -o ../globals.o

clean:
	rm -f $(OBJS) $(TARGET)
