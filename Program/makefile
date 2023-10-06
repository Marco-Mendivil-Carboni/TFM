NVCC = nvcc
NVCCFLAGS = -gencode arch=compute_35,code=sm_35 \
	-Wno-deprecated-gpu-targets -g
NVCCLINKS = -lcurand

bin/performsim: bin/util.o bin/chrdat.o bin/sugrid.o bin/chrsim.o \
		bin/performsim.o
	$(NVCC) $(NVCCFLAGS) $^ -o $@ $(NVCCLINKS)

bin/measuresim: bin/util.o bin/chrdat.o bin/sugrid.o bin/chrsim.o \
		bin/measuresim.o
	$(NVCC) $(NVCCFLAGS) $^ -o $@ $(NVCCLINKS)

bin/performsim.o: src/performsim.cu bin/chrsim.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/measuresim.o: src/measuresim.cu bin/chrsim.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/chrsim.o: src/chrsim.cu src/chrsim.cuh bin/chrdat.o bin/sugrid.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/sugrid.o: src/sugrid.cu src/sugrid.cuh bin/util.o src/vecops.cuh
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/chrdat.o: src/chrdat.cu src/chrdat.cuh bin/util.o src/vecops.cuh
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/util.o: src/util.cu src/util.cuh
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	rm bin/util.o bin/chrdat.o bin/sugrid.o bin/chrsim.o \
	bin/performsim.o bin/performsim bin/measuresim.o bin/measuresim
