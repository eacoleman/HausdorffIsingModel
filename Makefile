CC=g++
CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))
CC_FLAGS  := -fPIC -c -g -Wall `root-config --cflags` -02 -MMD
RT_FLAGS  := `root-config --cflags --glibs` -lMinuit -lMathMore -lMinuit2


all: runIsingModel 

test: libIsingModel.so
	$(CC) $(CC_FLAGS) scripts/testIsingModel.cpp -LlibIsingModel.so -o $@ $^ 

runIsingModel: libIsingModel.so
	$(CC) $(CC_FLAGS) scripts/runIsingModel.cpp -LlibIsingModel.so -o $@ $^ 

libIsingModel.so: $(OBJ_FILES) 
	$(CC) -fPIC -shared $(RT_FLAGS) -o $@ $^

obj/%.o: src/%.cpp
	g++ $(CC_FLAGS) -c -o $@ $<

clean:
	rm -f obj/*.o obj/*.d
	rm -f libIsingModel.so
