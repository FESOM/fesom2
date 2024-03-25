include Makefile.in

all:
	make metis
	make cleanfesom
	make fesom_ini
	make cleanfesom
	make fesom

metis:
	make -C ./lib/metis-5.1.0

fesom_ini:
	make -C ./src run_ini

fesom:
	make -C ./src run

cleanmetis:
	make -C ./lib/metis-5.1.0 clean

cleanfesom:
	make -C ./src clean

cleanall:
	make cleanmetis 
	make cleanfesom

clean:
	make cleanfesom
