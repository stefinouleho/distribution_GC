CFLAGS=-g -Wall
#CFLAGS=-O2 -Wall

#TEST= 28427
TEST= 15854
#TEST= 2361
#TEST= 16113
#TEST= 4672
#TEST= 28973
TAILLE=60
LC=7
#TEST= 28427

run: modelisation
	./modelisation  $(TAILLE)

val: modelisation
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./modelisation  ${TEST} $(TAILLE) $(LC)

modelisation: modelisation.o lecture_molecule_sdf.o fonctions_modelisation.o connexite.o
	gcc ${CFLAGS} modelisation.o lecture_molecule_sdf.o fonctions_modelisation.o connexite.o -o modelisation

modelisation.o: modelisation.c modelisation.h
	gcc ${CFLAGS} -c modelisation.c

fonctions_modelisation.o: fonctions_modelisation.c modelisation.h
		gcc ${CFLAGS} -c fonctions_modelisation.c

connexite.o: connexite.c modelisation.h
		gcc ${CFLAGS} -c connexite.c		
	
lecture_molecule_sdf.o: lecture_molecule_sdf.c lecture_molecule_sdf.h 
	gcc ${CFLAGS} -c lecture_molecule_sdf.c

clean: 
	rm -f modelisation
	rm -f modelisation.o
	rm -f lecture_molecule_sdf.o
	rm -f fonctions_modelisation.o
	rm -f connexite.o

