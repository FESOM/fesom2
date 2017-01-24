include Makefile.in

libmetis.a:
	make clean
	make -C ./GKlib
	make -C ./libmetis

clean: 
	rm -f ./libmetis/*.o  -f ./GKlib/*.o ./lib/libmetis.a

