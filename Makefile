all:
	g++ -s -o explicit_thr explicit_thr.cpp erand48.c _rand48.c  -march=native -funroll-loops -mfpmath=sse -msse2 -msse3 -mssse3 -ffast-math   -O3 -Wall -Wunused -Wno-format-y2k  -fno-exceptions -fno-strict-aliasing -fopenmp -mconsole  -luuid -lcomctl32 -static  -static-libgcc -static-libstdc++
	g++ -s -o explicit explicit.cpp erand48.c _rand48.c  -march=native -funroll-loops -mfpmath=sse -msse2 -msse3 -mssse3 -ffast-math   -O3 -Wall -Wunused -Wno-format-y2k  -fno-exceptions -fno-strict-aliasing -fopenmp -mconsole  -luuid -lcomctl32 -static  -static-libgcc -static-libstdc++
	g++ -s -o smallpt smallpt.cpp erand48.c _rand48.c  -march=native -funroll-loops -mfpmath=sse -msse2 -msse3 -mssse3 -ffast-math   -O3 -Wall -Wunused -Wno-format-y2k  -fno-exceptions -fno-strict-aliasing -fopenmp -mconsole  -luuid -lcomctl32 -static  -static-libgcc -static-libstdc++
	g++ -s -o ptreg ptreg.cpp erand48.c _rand48.c  -march=native -funroll-loops -mfpmath=sse -msse2 -msse3 -mssse3 -ffast-math   -O3 -Wall -Wunused -Wno-format-y2k  -fno-exceptions -fno-strict-aliasing -fopenmp -mconsole  -luuid -lcomctl32 -static  -static-libgcc -static-libstdc++
	
opt:
	g++ -s -o ptreg0-opt ptreg0-opt.cpp  -march=native -funroll-loops -mfpmath=sse -msse2 -msse3 -mssse3 -ffast-math   -O3 -Wall -Wunused -Wno-format-y2k  -fno-exceptions -fno-strict-aliasing -fopenmp -mconsole  -luuid -lcomctl32 -static  -static-libgcc -static-libstdc++
	
opt-new:
	g++ -s -o ptreg0-opt-new ptreg0-opt-new.cpp  -march=native -funroll-loops -mfpmath=sse -msse2 -msse3 -mssse3 -ffast-math   -O3 -Wall -Wunused -Wno-format-y2k  -fno-exceptions -fno-strict-aliasing -fopenmp -mconsole  -luuid -lcomctl32 -static  -static-libgcc -static-libstdc++
	
opt0:
	g++ -s -o ptreg0-opt0 ptreg0-opt0.cpp  -march=native -funroll-loops -mfpmath=sse -msse2 -msse3 -mssse3 -ffast-math   -O3 -Wall -Wunused -Wno-format-y2k  -fno-exceptions -fno-strict-aliasing -fopenmp -mconsole  -luuid -lcomctl32 -static  -static-libgcc -static-libstdc++
	
ex:
	g++ -s -o ptreg0-ex ptreg0-ex.cpp  -march=native -funroll-loops -mfpmath=sse -msse2 -msse3 -mssse3 -ffast-math   -O3 -Wall -Wunused -Wno-format-y2k  -fno-exceptions -fno-strict-aliasing -fopenmp -mconsole  -luuid -lcomctl32 -static  -static-libgcc -static-libstdc++
	
ptreg0:
	g++ -s -o ptreg0 ptreg0.cpp  -march=native -funroll-loops -mfpmath=sse -msse2 -msse3 -mssse3 -ffast-math   -O3 -Wall -Wunused -Wno-format-y2k  -fno-exceptions -fno-strict-aliasing -fopenmp -mconsole  -luuid -lcomctl32 -static  -static-libgcc -static-libstdc++
