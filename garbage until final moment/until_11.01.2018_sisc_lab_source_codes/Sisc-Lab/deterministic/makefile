output	:main.o multi.o	
			gcc -o output main.o multi.o -I /opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas -lpthread -lgfortran
multi.o		:multi.c
			gcc -c multi.c -I /opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas -lpthread -lgfortran

.PHONY:clean

clean:
		rm output main.o multi.o

