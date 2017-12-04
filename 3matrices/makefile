trial	:trial2.o multi.o	
			gcc -o trial trial2.o multi.o -I /opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas -lpthread -lgfortran
multi.o		:multi.c
			gcc -c multi.c -I /opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas -lpthread -lgfortran

.PHONY:clean

clean:
		rm trial trial2.o multi.o result.txt
