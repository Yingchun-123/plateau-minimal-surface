all: main 
Mesh2.o: Mesh2.cpp
	$(CXX) -c Mesh2.cpp

main: main.cpp Mesh2.o #GC.o
	$(CXX) main.cpp Mesh2.o -o main 


clean:
	-rm *.o *~ main 


