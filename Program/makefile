CC = g++
CCFLAGS = -g
CCLINKS = -lm

NVCC = nvcc
NVCCFLAGS = -gencode arch=compute_35,code=sm_35 -Wno-deprecated-gpu-targets -g
NVCCLINKS = -lcurand

bin/simulate: bin/main.o bin/chrsim.o bin/util.o
	$(NVCC) $(NVCCFLAGS) $^ -o $@ $(NVCCLINKS)

bin/main.o: src/main.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/chrsim.o: src/chrsim.cu src/chrsim.cuh
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/util.o: src/util.cpp src/util.hpp
	$(CC) $(CCFLAGS) -c $< -o $@

clean:
	rm bin/main.o bin/chrsim.o bin/util.o
