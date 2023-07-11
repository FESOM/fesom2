include Makefile.in

all:
	make metis
	make parms
	make cleanfesom
	make fesom_ini
	make cleanfesom
	make fesom

metis:
	make -C ./lib/metis-5.1.0

parms:
	make -C ./lib/parms

fesom_ini:
	make -C ./src run_ini

fesom:
	make -C ./src run

cleanparms:
	make -C ./lib/parms cleanall

cleanmetis:
	make -C ./lib/metis-5.1.0 clean

cleanfesom:
	make -C ./src clean

cleanall:
	make cleanparms 
	make cleanmetis 
	make cleanfesom

clean:
	make cleanfesom
