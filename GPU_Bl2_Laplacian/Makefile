#
# Makefile for Patus stencil benchmark
#
# Note: $(PATUS_*) variables will be automatically replaced by the
# required runtime files by Patus.
# find / -name cutil.h -type f -print 2>/dev/null
#

CC = nvcc
NVCCFLAGS = -O3 -arch=sm_13 -I/opt/NVIDIA_GPU_Computing_SDK/C/common/inc/ -L/opt/cluster/spool/cuda/cudatoolkit_3.0_linux_64_sled11.0/lib64
#NVCCFLAGS = -O0 -g -arch=sm_13 -I/home/christen/NVIDIA_GPU_Computing_SDK/C/common/inc

bench: kernel.cu driver.cu timer.c
	$(CC) $(NVCCFLAGS) -o $@ $+

clean:
	rm -rf *.o bench
