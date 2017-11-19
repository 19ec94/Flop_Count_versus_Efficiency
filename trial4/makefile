trial	:trial2.o allocate.o initialize.o initialize_zero.o multi.o	
			gcc -o trial trial2.o allocate.o initialize.o initialize_zero.o multi.o -I /opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas -lpthread -lgfortran
trial2.o	:trial2.c
			gcc -c trial2.c
allocate.o	:allocate.c
			gcc -c allocate.c 
initialize.o	:initialize.c
			gcc -c initialize.c
initialize_zero.o	:initialize_zero.c
			gcc -c initialize_zero.c
multi.o		:multi.c
			gcc -c multi.c -I /opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas -lpthread -lgfortran

.PHONY:clean

clean:
		rm trial trial2.o allocate.o initialize.o initialize_zero.o multi.o