#Variables

NVCC = nvcc
NVCCFLAGS = -gencode arch=compute_86,code=sm_86 -g

NVCCLINKS = -lcurand

#Targets

bin/makeallsim: bin/util.o bin/chrdat.o bin/makeallsim.o
	$(NVCC) $(NVCCFLAGS) $^ -o $@

bin/makeallsim.o: src/makeallsim.cu bin/performsim bin/analyzesim
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/performsim: bin/util.o bin/chrdat.o bin/sugrid.o bin/chrsim.o \
		bin/performsim.o
	$(NVCC) $(NVCCFLAGS) $^ -o $@ $(NVCCLINKS)

bin/measuresim: bin/util.o bin/chrdat.o bin/sugrid.o bin/chrsim.o \
		bin/measuresim.o
	$(NVCC) $(NVCCFLAGS) $^ -o $@ $(NVCCLINKS)

bin/analyzesim: bin/util.o bin/chrdat.o bin/chrana.o bin/analyzesim.o
	$(NVCC) $(NVCCFLAGS) $^ -o $@

bin/performsim.o: src/performsim.cu bin/chrsim.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/measuresim.o: src/measuresim.cu bin/chrsim.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/analyzesim.o: src/analyzesim.cu bin/chrana.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/chrsim.o: src/chrsim.cu src/chrsim.cuh bin/chrdat.o bin/sugrid.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/chrana.o: src/chrana.cu src/chrana.cuh bin/chrdat.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/chrdat.o: src/chrdat.cu src/chrdat.cuh bin/util.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/sugrid.o: src/sugrid.cu src/sugrid.cuh bin/util.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/util.o: src/util.cu src/util.cuh src/vect.cuh
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	rm -r bin/*
