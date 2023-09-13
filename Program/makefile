NVCC = nvcc
NVCCFLAGS = -gencode arch=compute_35,code=sm_35 \
	-Wno-deprecated-gpu-targets -g
NVCCLINKS = -lcurand

bin/simulate: bin/util.o bin/simulate.o \
	bin/chrdat.o bin/chrsim.o
	$(NVCC) $(NVCCFLAGS) $^ -o $@ $(NVCCLINKS)

bin/simulate.o: src/simulate.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/chrsim.o: src/chrsim.cu src/chrsim.cuh
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/chrdat.o: src/chrdat.cu src/chrdat.cuh
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/util.o: src/util.cu src/util.cuh
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	rm bin/util.o bin/simulate.o \
	bin/chrdat.o bin/chrsim.o
