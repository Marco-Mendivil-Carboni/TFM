#Variables

NVCC = nvcc
NVCCFLAGS = -gencode arch=compute_86,code=sm_86

NVCCLINKS = -lcurand

#Targets

all: bin/ccp-perform bin/ccp-profile bin/ccp-analyze

bin/ccp-perform: bin/util.o bin/chrdat.o bin/sugrid.o bin/chrsim.o \
		bin/ccp-perform.o
	$(NVCC) $(NVCCFLAGS) $^ -o $@ $(NVCCLINKS)

bin/ccp-profile: bin/util.o bin/chrdat.o bin/sugrid.o bin/chrsim.o \
		bin/ccp-profile.o
	$(NVCC) $(NVCCFLAGS) $^ -o $@ $(NVCCLINKS)

bin/ccp-analyze: bin/util.o bin/chrdat.o bin/sostat.o bin/chrana.o \
		bin/ccp-analyze.o
	$(NVCC) $(NVCCFLAGS) $^ -o $@

bin/ccp-perform.o: src/ccp-perform.cu bin/chrsim.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/ccp-profile.o: src/ccp-profile.cu bin/chrsim.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/ccp-analyze.o: src/ccp-analyze.cu bin/chrana.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/chrsim.o: src/chrsim.cu src/chrsim.cuh bin/chrdat.o bin/sugrid.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/chrana.o: src/chrana.cu src/chrana.cuh bin/chrdat.o bin/sostat.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/chrdat.o: src/chrdat.cu src/chrdat.cuh bin/util.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/sugrid.o: src/sugrid.cu src/sugrid.cuh bin/util.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/sostat.o: src/sostat.cu src/sostat.cuh bin/util.o
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

bin/util.o: src/util.cu src/util.cuh src/vect.cuh
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	rm -r bin/*
