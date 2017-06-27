OPT=-O3 -mfpmath=sse -ffast-math -fprefetch-loop-arrays -march=native -mtune=native

all: code docs clean

clean:
	rm -rf *.o *.mod *.aux *.bib *.log

code: global_env.o calendar.o cell.o newseed.o walker.o pool.o rw.o pooltest.o rw pooltest

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

pool.o: pool.F90
	gfortran ${OPT} -c pool.F90

rw.o: rw.F90
	gfortran ${OPT} -c rw.F90

pooltest.o: pooltest.F90
	gfortran ${OPT} -c pooltest.F90

rw: rw.o global_env.o calendar.o cell.o newseed.o walker.o
	gfortran ${OPT} -o rw rw.o global_env.o calendar.o cell.o newseed.o walker.o

pooltest: pool.o pooltest.o
	gfortran ${OPT} -o pooltest pooltest.o pool.o

docs: pdfbuilder
	./pdfbuilder
