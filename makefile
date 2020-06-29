src := Boris/*.cpp
 
obj := $(patsubst %.cpp,%.o,$(wildcard $(src)))
 
all: 
	g++ Boris/*.o -fopenmp -lsfml-graphics -lsfml-window -lsfml-system -lfftw3 -lX11 -o BorisLin
 
compile: $(obj)
 
%.o: %.cpp
	g++ -c -Ofast -std=c++17 -IBorisLib -IBorisCUDALib -fopenmp $< -o $@
 
clean: 
	rm -f $(obj) BorisLin
