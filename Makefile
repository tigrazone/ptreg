all:
	g++ -s -o explicit_thr explicit_thr.cpp erand48.c _rand48.c  -march=native -funroll-loops -mfpmath=sse -msse2 -msse3 -mssse3 -ffast-math   -O4 -Wall -Wunused -Wno-format-y2k  -fno-exceptions -fno-strict-aliasing -fopenmp -mconsole  -luuid -lcomctl32 -static  -static-libgcc -static-libstdc++
	g++ -s -o explicit explicit.cpp erand48.c _rand48.c  -march=native -funroll-loops -mfpmath=sse -msse2 -msse3 -mssse3 -ffast-math   -O4 -Wall -Wunused -Wno-format-y2k  -fno-exceptions -fno-strict-aliasing -fopenmp -mconsole  -luuid -lcomctl32 -static  -static-libgcc -static-libstdc++
	g++ -s -o smallpt smallpt.cpp erand48.c _rand48.c  -march=native -funroll-loops -mfpmath=sse -msse2 -msse3 -mssse3 -ffast-math   -O4 -Wall -Wunused -Wno-format-y2k  -fno-exceptions -fno-strict-aliasing -fopenmp -mconsole  -luuid -lcomctl32 -static  -static-libgcc -static-libstdc++
	g++ -s -o ptreg ptreg.cpp erand48.c _rand48.c  -march=native -funroll-loops -mfpmath=sse -msse2 -msse3 -mssse3 -ffast-math   -O4 -Wall -Wunused -Wno-format-y2k  -fno-exceptions -fno-strict-aliasing -fopenmp -mconsole  -luuid -lcomctl32 -static  -static-libgcc -static-libstdc++