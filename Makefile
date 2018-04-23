CC = gcc
SIM5LIB = sim5
RAY = Ray

CFLAGS = -I$(SIM5LIB) -L$(SIM5LIB) -Wall -O3 -w -fgnu89-inline 
LFLAGS = $(CFLAGS) -lm


default: sphere

%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

clean:
	rm -f *.o


sphere-src = \
    $(RAY)/$(SIM5LIB)/sim5lib2.c \
    $(RAY)/debugpro2.c \

sphere-obj = $(sphere-src:.c=.o)

sphere: $(sphere-obj)
	$(CC) $(sphere-obj) -o $@ $(LFLAGS) 


