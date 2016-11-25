CC=gcc
CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))
CC_FLAGS  := -Wall -g -std=gnu++11 -c `root-config --cflags --glibs` 
RT_FLAGS  := `root-config --cflags --glibs` -lMinuit -lMathMore -lMinuit2

all: test run

test: 
	$(CC) $(CC_FLAGS) src/testIsingModel.cpp -o testIsingModel 

run:
	$(CC) $(CC_FLAGS) src/runIsingModel.cpp -o testIsingModel 

clean:
	rm -f testIsingModel
	rm -f runIsingModel
