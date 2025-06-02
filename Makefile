all: sph

sph: sph.cpp
	g++ -o sph sph.cpp
#	g++ -O3 -march=native -ffast-math -fopt-info-vec-all -o sph sph.cpp 

clean:
	rm sph
