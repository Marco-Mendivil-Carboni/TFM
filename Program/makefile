CC = g++
CCFLAGS = -g
CCLINKS = -lm

NVCC = nvcc
NVCCFLAGS = -gencode arch=compute_35,code=sm_35 -Wno-deprecated-gpu-targets -g
NVCCLINKS = -lcurand

bin/simulate: bin/simulate.o bin/chromatin.o bin/utilities.o
	$(NVCC) $(NVCCFLAGS) $^ -o $@ $(NVCCLINKS)

bin/simulate.o: src/simulate.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/chromatin.o: src/chromatin.cu src/chromatin.cuh
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/utilities.o: src/utilities.cpp src/utilities.hpp
	$(CC) $(CCFLAGS) -c $< -o $@

clean:
	rm bin/simulate.o bin/chromatin.o bin/utilities.o
