all:
	g++ -s -o ptreg ptreg.cpp erand48.c _rand48.c  -march=native -funroll-loops -mfpmath=sse -msse2 -msse3 -mssse3 -ffast-math   -O4 -Wall -Wunused -Wno-format-y2k  -fno-exceptions -fno-strict-aliasing -fopenmp -mconsole  -static  -static-libgcc -static-libstdc++
