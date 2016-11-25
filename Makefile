CC=c++
CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))
CC_FLAGS  := -fPIC -std=c++11 -c -g -Wall `root-config --cflags --glibs` -MMD 
RT_FLAGS  := `root-config --cflags --glibs` -lMinuit -lMathMore -lMinuit2


all: runIsingModel 

test: libIsingModel.so
	$(CC) $(CC_FLAGS) -I$(ROOTSYS)/include -L$(ROOTSYS)/lib scripts/testIsingModel.cpp -LlibIsingModel.so -o $@ $^ 

runIsingModel: libIsingModel.so
	$(CC) $(CC_FLAGS) scripts/runIsingModel.cpp -LlibIsingModel.so -o $@ $^ 

libIsingModel.so: $(OBJ_FILES) 
	$(CC) -fPIC -shared $(RT_FLAGS) -o $@ $^

obj/%.o: src/%.cpp
	$(CC) $(CC_FLAGS) -o $@ $<

clean:
	rm -f obj/*.o obj/*.d
	rm -f libIsingModel.so
