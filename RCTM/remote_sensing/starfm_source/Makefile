#############################################################################
# makefile for StarFM
#
# makefile name: Makefile (for linux)
#
# Revision 1.2.1   10/2014  enabled OpenMP for multi-thread processing
# Revision 1.1     01/2008  initial makefile
#
##!END
#############################################################################
#CC = gcc -Wall -O2
#CC = gcc -Wall -O2 -fopenmp -DDEBUG
CC = gcc -Wall -O2 -fopenmp

TARGET = StarFM.exe

OBJ = StarFM_main.o StarFM_compute.o StarFM_util.o StarFM_alloc.o

all: $(TARGET)

$(OBJ) :StarFM.h

# Make the process
$(TARGET) : $(OBJ)
	$(CC) $(OBJ) -lm  -o $(TARGET)

clean:
	rm -f *.o

#******************* End of make file *******************************


