MPICC=mpic++
RAND_FLAG = -I ./Random123/include/ -O2

pattern : CPTF_pattern.cpp
	$(MPICC) CPTF_pattern.cpp -o $@

seedfind : CPTF_seedfind.cpp
	$(MPICC) CPTF_seedfind.cpp -o $@
