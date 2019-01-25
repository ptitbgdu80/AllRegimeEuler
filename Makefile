# Compilateur utilisé
CC=g++

# Les options de compilation
FLAGS = -std=c++11 -Wall

# Le nom de l'exécutable
PROG = test
PROG1D = 1D
PROG2DPC = 2DPC
PROGR2D = R2D

# Les fichiers source à compiler
SRC = main.cc AREuler.cpp AREuler2D.cpp
SRC1D0 = 1Dtheta0.cc AREuler.cpp AREuler2D.cpp
SRC1D1 = 1Dtheta1.cc AREuler.cpp AREuler2D.cpp
SRC1D2 = 1Dtheta2.cc AREuler.cpp AREuler2D.cpp
SRC2DPC0 = 2Dtheta0PC.cc AREuler.cpp AREuler2D.cpp
SRC2DPC1 = 2Dtheta1PC.cc AREuler.cpp AREuler2D.cpp
SRC2DPC2 = 2Dtheta2PC.cc AREuler.cpp AREuler2D.cpp
SRCR2D0 = Riemann2Dtheta0.cc AREuler.cpp AREuler2D.cpp
SRCR2D1 = Riemann2Dtheta1.cc AREuler.cpp AREuler2D.cpp
SRCR2D2 = Riemann2Dtheta2.cc AREuler.cpp AREuler2D.cpp
SRCR2D3 = Riemann2DRusa.cc AREuler.cpp AREuler2D.cpp

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

# La commande complète : compile seulement si un fichier a été modifié
$(PROG2DPC) : $(SRC2DPC0) $(SRC2DPC1) $(SRC2DPC2)
	$(CC) $(FLAGS) $(SRC2DPC0) -o run0
	./run0 &
	$(CC) $(FLAGS) $(SRC2DPC1) -o run1
	./run1 &
	$(CC) $(FLAGS) $(SRC2DPC2) -o run2
	./run2 &

# La commande complète : compile seulement si un fichier a été modifié
$(PROGR2D) : $(SRCR2D0) $(SRCR2D1) $(SRCR2D2) $(SRCR2D3)
	$(CC) $(FLAGS) $(SRCR2D0) -o run0
	./run0 &
	$(CC) $(FLAGS) $(SRCR2D1) -o run1
	./run1 &
	$(CC) $(FLAGS) $(SRCR2D2) -o run2
	./run2 &
	$(CC) $(FLAGS) $(SRCR2D3) -o run3
	./run3 &

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ $(PROG1D)
