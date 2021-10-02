#if cuda architecture not specified then assume 50, i.e. -arch=sm_50
#arch value for required architecture should be:
#arch=50 for Maxwell
#arch=60 for Pascal
#arch=70 for Volta (and Turing)
#arch=80 for Ampere
#example: $ make configure arch=70
ifndef arch
	arch = $(file < arch.txt)
endif

#if you want to compile cuda code in double precision then $ make configure sprec=0
ifndef sprec
	sprec = 1
endif

#Boris program version
BVERSION := 340

#Working directories
BORIS_DATA_DIR := Boris_Data
BORIS_SIM_DIR := $(BORIS_DATA_DIR)/Simulations

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
	rm -f $(OBJ_FILES) $(CUOBJ_FILES) BorisLin

#compile only cpp files
cpp: $(OBJ_FILES)
	@echo Done
	
#compile only cu files
cuda: $(CUOBJ_FILES)
	@echo Done
	
#configure CUDA compilation first: architecture and float precision
configure:
	$(file > BorisCUDALib/cuBLib_Flags.h,#pragma once)
	$(file >> BorisCUDALib/cuBLib_Flags.h,)
	$(file >> BorisCUDALib/cuBLib_Flags.h,)
	$(file >> BorisCUDALib/cuBLib_Flags.h,#define __CUDA_ARCH__ $(arch)0)
	$(file >> BorisCUDALib/cuBLib_Flags.h,)
	$(file >> BorisCUDALib/cuBLib_Flags.h,)
	$(file >> BorisCUDALib/cuBLib_Flags.h,#define SINGLEPRECISION $(sprec))
	$(file > arch.txt,$(arch))
	mkdir -p $(OBJ_DIR)
	mkdir -p $(CUOBJ_DIR)
	@echo Configured for -arch=sm_$(arch) and SINGLEPRECISION = $(sprec)
	
#compile both cpp and cu files
compile: $(OBJ_FILES) $(CUOBJ_FILES)
	@echo Done
  
install:
	nvcc -arch=sm_$(arch) -dlink -w $(CUOBJ_DIR)/*.o -o $(CUOBJ_DIR)/rdc_link.o 
	g++ $(OBJ_DIR)/*.o $(CUOBJ_DIR)/*.o -fopenmp -ltbb -lsfml-graphics -lsfml-window -lsfml-system -lfftw3 -lX11 -lcudart -lcufft -lcudadevrt -o BorisLin
	rm -f $(OBJ_FILES) $(CUOBJ_FILES) $(CUOBJ_DIR)/rdc_link.o
	mkdir -p ~/Documents/$(BORIS_DATA_DIR)
	mkdir -p ~/Documents/$(BORIS_SIM_DIR)
	cp -r BorisPythonScripts/ ~/Documents/$(BORIS_DATA_DIR)
	cp -r Examples/ ~/Documents/$(BORIS_DATA_DIR)
	cp -f BorisMDB.txt ~/Documents/$(BORIS_DATA_DIR)/BorisMDB.txt
	cp -f Manual/BorisManual-v$(BVERSION).pdf ~/Documents/$(BORIS_DATA_DIR)/BorisManual-v$(BVERSION).pdf
	@echo Done
 
Boris/Boris_o/%.o: Boris/%.cpp
	g++ -c -Ofast -std=c++17 -IBorisLib -IBorisCUDALib -fopenmp $< -o $@

Boris/Boris_cuo/%.o: Boris/%.cu
	nvcc -rdc=true -c -std=c++14 -IBorisLib -IBorisCUDALib -w -arch=sm_$(arch) $< -o $@
	
