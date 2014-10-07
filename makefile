hidro0: hidro-mcv.o kernel.o SPH-state.o SPH-lista.o SPH-ODEsolver.o utilities.o src/hydro/SPHtypes.h src/hydro/SPH-lista.h src/hydro/SPH-ODEsolver.h src/hydro/kernel.h
	gcc -o hidro0 obj/hidro-mcv.o obj/utilities.o obj/kernel.o obj/SPH-state.o obj/SPH-lista.o obj/SPH-ODEsolver.o -O3 -lm

hidro1: hidro-mcv.o kernel.o SPH-state.o SPH-lista.o SPH-ODEsolver.o utilities.o src/hydro/SPHtypes.h src/hydro/SPH-lista.h src/hydro/SPH-ODEsolver.h src/hydro/kernel.h
	gcc -o hidro1 obj/hidro-mcv.o obj/utilities.o obj/kernel.o obj/SPH-state.o obj/SPH-lista.o obj/SPH-ODEsolver.o -O3 -lm

hidro2: hidro-mcv.o kernel.o SPH-state.o SPH-lista.o SPH-ODEsolver.o utilities.o src/hydro/SPHtypes.h src/hydro/SPH-lista.h src/hydro/SPH-ODEsolver.h src/hydro/kernel.h
	gcc -o hidro2 obj/hidro-mcv.o obj/utilities.o obj/kernel.o obj/SPH-state.o obj/SPH-lista.o obj/SPH-ODEsolver.o -O3 -lm

hidro3: hidro-mcv.o kernel.o SPH-state.o SPH-lista.o SPH-ODEsolver.o utilities.o src/hydro/SPHtypes.h src/hydro/SPH-lista.h src/hydro/SPH-ODEsolver.h src/hydro/kernel.h
	gcc -o hidro3 obj/hidro-mcv.o obj/utilities.o obj/kernel.o obj/SPH-state.o obj/SPH-lista.o obj/SPH-ODEsolver.o -O3 -lm	

print:  print_extern.o kernel.o SPH-state.o SPH-lista.o SPH-ODEsolver.o utilities.o src/hydro/SPHtypes.h src/hydro/SPH-lista.h src/hydro/SPH-ODEsolver.h src/hydro/kernel.h
	gcc -o print obj/print_extern.o obj/utilities.o obj/kernel.o obj/SPH-state.o obj/SPH-lista.o obj/SPH-ODEsolver.o -O3 -lm

print_extern.o: src/hydro/print_extern.c src/hydro/SPH-lista.h src/hydro/SPH-state.h src/hydro/kernel.h src/hydro/SPHtypes.h src/hydro/SPH-ODEsolver.h src/hydro/utilities.h
	gcc -c src/hydro/print_extern.c -O3 -ansi -o obj/print_extern.o

kernel.o: src/hydro/kernel.c
	gcc -c src/hydro/kernel.c -O3 -ansi -o obj/kernel.o

SPH-state.o: src/hydro/SPH-state.c src/hydro/SPHtypes.h
	gcc -c src/hydro/SPH-state.c -O3 -ansi -o obj/SPH-state.o

SPH-lista.o: src/hydro/SPH-lista.c src/hydro/SPHtypes.h src/hydro/SPH-state.h
	gcc -c src/hydro/SPH-lista.c -O3 -ansi -o obj/SPH-lista.o

SPH-ODEsolver.o: src/hydro/SPH-ODEsolver.c src/hydro/SPHtypes.h
	gcc -c src/hydro/SPH-ODEsolver.c -O3 -ansi -o obj/SPH-ODEsolver.o

utilities.o: src/hydro/utilities.c src/hydro/SPH-lista.h src/hydro/SPH-state.h src/hydro/kernel.h src/hydro/SPHtypes.h src/hydro/SPH-ODEsolver.h src/hydro/utilities.h
	gcc -c src/hydro/utilities.c -O3 -ansi -o obj/utilities.o

hidro-mcv.o: src/hydro/hidro-mcv.c src/hydro/SPH-lista.h src/hydro/SPH-state.h src/hydro/kernel.h src/hydro/SPHtypes.h src/hydro/SPH-ODEsolver.h src/hydro/utilities.h
	gcc -c src/hydro/hidro-mcv.c -O3 -ansi -o obj/hidro-mcv.o
