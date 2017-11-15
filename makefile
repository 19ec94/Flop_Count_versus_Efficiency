trial	:trial2.o allocate.o initialize.o allocate_zero.o
			gcc -o trial trial2.o allocate.o initialize.o allocate_zero.o /home/edilbert/library/OpenBLAS-0.2.20/libopenblas.a -lpthread -lgfortran

trial2.o:trial2.c
			gcc -c trial2.c 
initialize.o:initialize.c
			gcc -c initialize.c
allocate.o:allocate.c
			gcc -c allocate.c
allocate_zero.o:allocate_zero.c
			gcc -c allocate_zero.c

.PHONY:clean

clean:
		rm trial trial2.o initialize.o allocate.o allocate_zero.o