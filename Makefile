all: sph

sph: sph.cpp
	g++ -o sph sph.cpp

clean:
	rm sph