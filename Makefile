OBJS  = dis.o        \
	var.o        \
	dif.o        \
	evo.o        \
	grd.o        \
	ini.o        \
	slp.o        \
	str.o        \
	ino.o        \
	alo.o        \
	gsf.o        \
	lin.o        \
	fft.o

FC=gfortran

FFLAGS= -O3 -I"/usr/local/include" -L"/usr/local/lib"

DIS:	${OBJS}  Makefile
	${FC} -o $@ ${OBJS} ${FFLAGS} -lfftw3_omp -lfftw3 -lm -fopenmp 

var.mod: var.o var.f
	 ${FC} -c var.f 
var.o:   var.f
	 ${FC} -c var.f 

dis.o:   var.mod dis.f
	 ${FC} -c dis.f -fopenmp -ffree-form 

dif.o:   var.mod dif.f
	 ${FC} -c dif.f -fopenmp  
evo.o:   var.mod evo.f
	 ${FC} -c evo.f -fopenmp  
grd.o:   var.mod grd.f
	 ${FC} -c grd.f -fopenmp
ini.o:   var.mod ini.f
	 ${FC} -c ini.f -fopenmp -ffree-form
slp.o:   var.mod slp.f
	 ${FC} -c slp.f
str.o:   var.mod str.f
	 ${FC} -c str.f -fopenmp 
ino.o:   var.mod ino.f
	 ${FC} -c ino.f -fopenmp 
alo.o:   var.mod alo.f
	 ${FC} -c alo.f -ffree-form
gsf.o:   var.mod gsf.f
	 ${FC} -c gsf.f -fopenmp
lin.o:   var.mod lin.f
	 ${FC} -c lin.f 
fft.o:   fft.f
	 ${FC} -c fft.f ${FFLAGS} -lfftw3_omp -lfftw3 -lm -fopenmp -ffree-form 
	 
clean:
	rm -f *.o *~ *.mod 

cleanall:
	rm -f *.o *~ DIS *.dat fort.* *.prf *.med *.nh *.dis *.out *.mod
