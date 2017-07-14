OPT=-O3 -mfpmath=sse -ffast-math -fprefetch-loop-arrays -march=native -mtune=native

all: code clean

clean:
	rm -rf *.o *.mod *.aux *.bib *.log

code: global_env.o calendar.o cell.o newseed.o walker.o blockrw.o blockrw

global_env.o: global_env.F90
	gfortran ${OPT} -c global_env.F90

calendar.o: calendar.F90
	gfortran ${OPT} -c calendar.F90

cell.o: cell.F90
	gfortran ${OPT} -c cell.F90

newseed.o: newseed.F90
	gfortran ${OPT} -c newseed.F90

walker.o: walker.F90
	gfortran ${OPT} -c walker.F90

blockrw.o: blockrw.F90
	gfortran ${OPT} -c blockrw.F90

blockrw: blockrw.o global_env.o calendar.o cell.o newseed.o walker.o
	gfortran ${OPT} -o blockrw blockrw.o global_env.o calendar.o cell.o newseed.o walker.o

