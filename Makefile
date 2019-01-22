# Compilateur utilisé
CC=g++

# Les options de compilation
FLAGS = -std=c++11 -Wall

# Le nom de l'exécutable
PROG = test
PROG1D = 1D

# Les fichiers source à compiler
SRC = main.cc AREuler.cpp AREuler2D.cpp
SRC1D0 = 1Dtheta0.cc AREuler.cpp AREuler2D.cpp
SRC1D1 = 1Dtheta1.cc AREuler.cpp AREuler2D.cpp
SRC1D2 = 1Dtheta2.cc AREuler.cpp AREuler2D.cpp

# La commande complète : compile seulement si un fichier a été modifié
$(PROG) : $(SRC)
	$(CC) $(FLAGS) $(SRC) -o run


# La commande complète : compile seulement si un fichier a été modifié
$(PROG1D) : $(SRC1D0) $(SRC1D1) $(SRC1D2)
	$(CC) $(FLAGS) $(SRC1D0) -o run0
	./run0
	$(CC) $(FLAGS) $(SRC1D1) -o run1
	./run1
	$(CC) $(FLAGS) $(SRC1D2) -o run2
	./run2
	rm -f *.o *~ $(PROG1D)
	# gnuplot finalSolution1D.gp


# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ $(PROG1D)
