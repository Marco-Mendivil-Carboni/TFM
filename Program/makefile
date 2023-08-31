CC = g++
CCFLAGS = -g
CCLINKS = -lm

NVCC = nvcc
NVCCFLAGS = -gencode arch=compute_35,code=sm_35 -Wno-deprecated-gpu-targets -g
NVCCLINKS = -lcurand

bin/simulate: bin/main.o bin/chrsim.o bin/util.o
	$(NVCC) $(NVCCFLAGS) $^ -o $@ $(NVCCLINKS)

bin/main.o: main.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/chrsim.o: src/chrsim.cu inc/chrsim.cuh
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/util.o: src/util.cpp inc/util.hpp
	$(CC) $(CCFLAGS) -c $< -o $@

clean:
	rm bin/main.o bin/chrsim.o bin/util.o
