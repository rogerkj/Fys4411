# Standard makefile

CC= c++ -Wall
PROG = main
LIB1 = lib
LIB2 = mcintegrator

LINK =  -llapack -lblas -larmadillo

${PROG} : 	${PROG}.o ${LIB1}.o  ${LIB2}.o 
		${CC} ${PROG}.o ${LIB1}.o ${LIB2}.o -o ${PROG} ${LINK}


clean:
	rm -f ${PROG} ${PROG}.o
	rm -f ${LIB} ${LIB}.o