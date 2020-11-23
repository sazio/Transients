CC = gcc 
CFLAGS= -O3 -ffast-math
LDFLAGS=-lm -lcvode -lcvodeserial

OBJECTS =  hcrt_cvode.o LC_SFO.o rgenerator.o

LC: $(OBJECTS)
	$(CC) $(CFLAGS) -o lc  LC_SFO.o rgenerator.o $(LDFLAGS)

HCRT: $(OBJECTS)
	$(CC) $(CFLAGS) -o hcrt hcrt_cvode.o rgenerator.o $(LDFLAGS)

all : LC HCRT

clean :
	rm -f *.o LC HCRT
