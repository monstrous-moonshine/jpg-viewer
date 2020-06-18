objects = main.o huffman.o image.o jpeg.o

app : $(objects) -lsfml-graphics -lsfml-window -lsfml-system
	g++ -o $@ $^

huffman.o : huffman.h
image.o   : image.h
jpeg.o    : jpeg.h

.PHONY : clean
clean :
	-rm app $(objects)
