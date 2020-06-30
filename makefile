src_cpu := Boris/*.cpp
src_gpu := Boris/*.cu
 
cpu := $(patsubst %.cpp,%.o,$(wildcard $(src_cpu)))
gpu := $(patsubst %.cu,%.co,$(wildcard $(src_gpu)))
 
.PHONY: clean
clean: 
	rm -f $(cpu) $(gpu) BorisLin
 
install: 
	g++ Boris/*.o Boris/*.co -fopenmp -lsfml-graphics -lsfml-window -lsfml-system -lfftw3 -lX11 -lcudart -lcufft -o BorisLin
	rm -f $(cpu) $(gpu)
 
compile: $(cpu) $(gpu)

cuda: $(gpu)
 
%.o: %.cpp
	g++ -c -Ofast -std=c++17 -IBorisLib -IBorisCUDALib -fopenmp $< -o $@

%.co: %.cu
	nvcc -c -std=c++14 -IBorisLib -IBorisCUDALib -w -arch=sm_50 $< -o $@
