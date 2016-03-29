#CC = g++ -O0 #valgrind
CC = g++ -O3 
DEBUG = 0
VERBOSE = 0
GPROF = 0

#INCDIRS = -I. -I${HOME}/boost_1_35_0
INCDIRS = -I. 
LIBS = -lm

ifeq (${DEBUG}, 1)
	GDB = -ggdb 
else
	GDB =
endif

ifeq (${GPROF}, 1)
	GPROF_PRM = -pg -g 
else
	GPROF_PRM =
endif

#-O3		For optimization
#-Wall		To turn on "all warnings"
#-m32		The 32-bit environment sets int, long and pointer to 32 bits. 
#
CFLAGS = -Wall -m32 ${GDB} ${GPROF_PRM} -D DEBUG=${DEBUG} -D VERBOSE=${VERBOSE} ${INCDIRS}

OBJS_PALIGN  = palign_main.cpp nwalign.o Input.o DisplayResults.o dataset.o symbols.o mt19937ar.o

all: palign 

.c.o .cpp.o: 
	${CC} ${CFLAGS} -c $<


palign: ${OBJS_PALIGN}
	${CC} ${CFLAGS} -o palign.out ${OBJS_PALIGN} ${LIBS}
	mv palign.out ../

clean: 
	@ \rm -f *.o depend

depend: 
	${CC} ${CFLAGS} -MM *.c *.cpp > depend

include depend
