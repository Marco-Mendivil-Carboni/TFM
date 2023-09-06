CC = g++
CCFLAGS = -g
CCLINKS = -lm

NVCC = nvcc
NVCCFLAGS = -gencode arch=compute_35,code=sm_35 -Wno-deprecated-gpu-targets -g
NVCCLINKS = -lcurand

# EXT_INC_DIR = -I/usr/include/tirpc
# EXT_LIB_DIR = 
# EXT_LINKS = -ltirpc

bin/simulate: bin/simulate.o bin/chrsim.o bin/util.o
	$(NVCC) $(NVCCFLAGS) $^ -o $@ $(NVCCLINKS)
#	 $(NVCC) $(NVCCFLAGS) $^ -o $@ $(NVCCLINKS) $(EXT_LINKS)

bin/simulate.o: src/simulate.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/chrsim.o: src/chrsim.cu src/chrsim.cuh
	$(NVCC) $(NVCCFLAGS) -c $< -o $@
#	 $(NVCC) $(NVCCFLAGS) $(EXT_INC_DIR) -c $< -o $@

bin/util.o: src/util.cpp src/util.hpp
	$(CC) $(CCFLAGS) -c $< -o $@

clean:
	rm bin/simulate.o bin/chrsim.o bin/util.o
