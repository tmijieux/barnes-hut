TARGET=particle
CXXFLAGS=-std=c++11 -g -fopenmp -Wall -Wextra # $(shell pkg-config --cflags glib-2.0) -fdiagnostics-color=auto
LIBS=-lm # -$(shell pkg-config --libs glib-2.0) 
LDFLAGS=-fopenmp
GENGETOPT=gengetopt
CXX=mpic++

ifdef DEBUG
CXXFLAGS+=-ggdb -O0 -DDEBUG=1
else
CXXFLAGS+=-O3
endif

OBJ=	main.o \
	../perf/perf.o \
	site.o \
	particle.o \
	proc.o

DEP=$(OBJ:.o=.d)

all: $(TARGET)

-include $(DEP)


particle: $(OBJ)
	$(CXX) $^ -o $@ $(LDFLAGS) $(LIBS)

%.o: %.cpp
	@$(CXX) -MM $(CXXFLAGS) $*.cpp > $*.d
	$(CXX) -c $(CXXFLAGS) $*.cpp -o $*.o

%.o: %.c
	@$(CXX) -MM $(CXXFLAGS) $*.c > $*.d
	$(CXX) -c $(CXXFLAGS) $*.c -o $*.o

clean:
	$(RM) $(OBJ) $(DEP) *.d *.o *~

mrproper: clean
	$(RM) $(TARGET)

genopt: particle.ggo
	$(GENGETOPT) -u"INPUT FILES" < particle.ggo

