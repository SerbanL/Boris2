#if cuda architecture not specified then assume 50, i.e. -arch=sm_50
#arch value for required architecture should be:
#arch=50 for Maxwell
#arch=60 for Pascal
#arch=70 for Volta (and Turing)
#arch=80 for Ampere
#arch=90 for Ada (and Hopper)
#example: $ make configure arch=80 sprec=1 python=3.8 cuda=12.0
ifndef arch
	arch = $(file < arch.txt)
endif

ifndef python
	python = $(file < python.txt)
endif

ifndef cuda
	cuda = $(file < cuda.txt)
endif

#if you want to compile cuda code in double precision then $ make configure sprec=0
ifndef sprec
	sprec = 1
endif

#Boris program version
BVERSION := 380

#Source directories for make
SRC_DIR := Boris
OBJ_DIR := Boris/Boris_o
CUOBJ_DIR := Boris/Boris_cuo
SRC_CPP_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_CPP_FILES))
SRC_CU_FILES := $(wildcard $(SRC_DIR)/*.cu)
CUOBJ_FILES := $(patsubst $(SRC_DIR)/%.cu,$(CUOBJ_DIR)/%.o,$(SRC_CU_FILES))

.PHONY: clean
clean: 
	rm -f $(OBJ_FILES) $(CUOBJ_FILES) $(CUOBJ_DIR)/rdc_link.o BorisLin

#compile only cpp files
cpp: $(OBJ_FILES)
	@echo Done
	
#compile only cu files
cuda: $(CUOBJ_FILES)
	@echo Done
	
#configure CUDA compilation first: architecture and float precision. Also set Python version and CUDA Toolkit version.
configure:
	$(file > BorisCUDALib/cuBLib_Flags.h,#pragma once)
	$(file >> BorisCUDALib/cuBLib_Flags.h,)
	$(file >> BorisCUDALib/cuBLib_Flags.h,)
	$(file >> BorisCUDALib/cuBLib_Flags.h,#define __CUDA_ARCH__ $(arch)0)
	$(file >> BorisCUDALib/cuBLib_Flags.h,)
	$(file >> BorisCUDALib/cuBLib_Flags.h,)
	$(file >> BorisCUDALib/cuBLib_Flags.h,#define SINGLEPRECISION $(sprec))
	$(file > arch.txt,$(arch))
	$(file > python.txt,$(python))
	$(file > cuda.txt,$(cuda))
	mkdir -p $(OBJ_DIR)
	mkdir -p $(CUOBJ_DIR)
	@echo Configured for -arch=sm_$(arch) and SINGLEPRECISION = $(sprec). Python version $(python). CUDA Toolkit version $(cuda).
	
#compile both cpp and cu files
compile: $(OBJ_FILES) $(CUOBJ_FILES)
	@echo Done
  
install:
	nvcc -arch=sm_$(arch) -dlink -w $(CUOBJ_DIR)/*.o -o $(CUOBJ_DIR)/rdc_link.o
	g++ $(OBJ_DIR)/*.o $(CUOBJ_DIR)/*.o -fopenmp -lpython$(python) -ltbb -lfftw3 -lX11 -L/usr/local/cuda-$(cuda)/targets/x86_64-linux/lib/ -lcudart -lcufft -lcudadevrt -o BorisLin
	#rm -f $(OBJ_FILES) $(CUOBJ_FILES) $(CUOBJ_DIR)/rdc_link.o
	@echo Done
 
#for python3.8 make sure to get dev version : sudo apt-get install python3.8-dev
Boris/Boris_o/%.o: Boris/%.cpp
	g++ -I/usr/local/cuda-$(cuda)/targets/x86_64-linux/include/ -c -Ofast -std=c++17 -I/usr/include/python$(python)/ -IBorisLib -IBorisCUDALib -fopenmp $< -o $@

Boris/Boris_cuo/%.o: Boris/%.cu
	nvcc -I/usr/local/cuda-$(cuda)/targets/x86_64-linux/include/ -rdc=true -c -std=c++14 -IBorisLib -IBorisCUDALib -w -arch=sm_$(arch) $< -o $@