src := Boris/*.cpp
 
obj := $(patsubst %.cpp,%.o,$(wildcard $(src)))
 
all: 
	g++ Boris/*.o -fopenmp -lsfml-graphics -lsfml-window -lsfml-system -lfftw3 -lX11 -o Boris
 
compile: $(obj)
 
%.o: %.cpp
	g++ -c -O0 -std=c++17 -IBorisLib -IBorisCUDALib -fopenmp $< -o $@
 
clean: 
	rm -f $(obj) Boris
