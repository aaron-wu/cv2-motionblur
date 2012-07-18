export GXX = g++ -O3 -Wall
export INCLUDE = .
export LIBRARY = .
export HEADERS = $(notdir $(wildcard ${INCLUDE}/**/*))
SOURCES = $(notdir $(wildcard *.cpp))
EXCLUDES = $(patsubst %.cpp, %.o, $(notdir main.cpp))
export OBJECTS = $(filter-out ${EXCLUDES}, $(patsubst %.cpp, %.o, $(SOURCES)))
PROG = motionblur
INCLUDE_ARGS = ${INCLUDE:%=-I%}
LIBRARY_ARGS = ${LIBRARY:%=-L%}
MAIN = main.cpp

all: prog

clean:
	rm -f *.o
	rm -f *.a
	rm -f ${PROG}

%.o: ${SRC}/%.cpp	
	$(GXX) -c $< ${INCLUDE_ARGS}

prog: ${OBJECTS}
	$(GXX) -o ${PROG} ${INCLUDE_ARGS} ${LIBRARY_ARGS} $^ ${MAIN}
	
compile: ${OBJECTS}

lib: compile
	ar -rv libbfr.a *.o

.PRECIOUS: %.o %.a

.PHONY: clean
