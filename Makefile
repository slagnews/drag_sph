CXX = g++
CXXFLAGS = -std=c++17 -Wall -O3 -fopenmp -Iglad/include -Iinclude
LIBS = -lglfw -ldl -lGL

HEADERS = include/sph.h include/shaders.h

TARGET = sph

SRC = src/main.cpp src/shaders.cpp glad/src/glad.c src/sph.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LIBS)

clean:
	rm -f $(TARGET)
