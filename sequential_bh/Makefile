TARGET=particle
CFLAGS=-std=gnu99 -g -fopenmp -Wall -Wextra # $(shell pkg-config --cflags glib-2.0) -fdiagnostics-color=auto
LIBS=-lm # -$(shell pkg-config --libs glib-2.0) 
LDFLAGS=-fopenmp
GENGETOPT=gengetopt
CC=gcc

ifdef DEBUG
CFLAGS+=-ggdb -O0 -DDEBUG=1
else
CFLAGS+=-O3
endif

SRC=	main.c \
	../perf/perf.c \
	site.c \
	particle.c

OBJ=$(SRC:.c=.o)
DEP=$(SRC:.c=.d)

all: $(TARGET)

-include $(DEP)


particle: $(OBJ)
	$(CC) $^ -o $@ $(LDFLAGS) $(LIBS)

%.o: %.c
	@$(CC) -MM $(CFLAGS) $*.c > $*.d
	$(CC) -c $(CFLAGS) $*.c -o $*.o

clean:
	$(RM) $(OBJ) $(DEP) *.d *.o *~

mrproper: clean
	$(RM) $(TARGET)

genopt: particle.ggo
	$(GENGETOPT) -u"INPUT FILES" < particle.ggo

