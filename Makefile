CC = gcc
RM = rm -rf
GSL_DIR = /opt/gsl/gcc-7.5.0/2.6/
CFLAGS = -g -O2 -Wall -mcmodel=medium -I$(GSL_DIR)/include -Wno-unused-variable
LDFLAGS = -L$(GSL_DIR)/lib 
LINKS = -lz -lm -lgsl -lgslcblas
OBJS = gf.o
TARGETS = gf

.PHONY: all clean dep
.SUFFIXES : .c .o

.c .o :
	$(CC) $(CFLAGS) -c $<

all : $(TARGETS)

gf  : gf.o
	$(CC) $(LDLIBS) $(LDFLAGS) -o $@ gf.o $(LINKS)

clean :
	$(RM) *.o
	$(RM) $(TARGET)

dep :
	$(CC) $(CFLAGS) -M $(OBJS:.o=.c) 

