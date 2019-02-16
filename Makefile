CC:=gcc
CFLAGS:=-Wall -O3 -msse4 -g
OBJ_DIR:=obj
BIN_DIR:=bin
SRC_DIR:=src
LIB_DIR:=lib

OBJ=gssw.o

.PHONY:all clean cleanlocal test

all: $(LIB_DIR)/libgssw.a

$(OBJ_DIR)/$(OBJ):$(SRC_DIR)/gssw.h $(SRC_DIR)/gssw.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c -o $@ $(SRC_DIR)/gssw.c

$(LIB_DIR)/libgssw.a:$(OBJ_DIR)/$(OBJ)
	@mkdir -p $(@D)
	ar rvs $@ $<
	
cleanlocal:
	$(RM) -r lib/
	$(RM) -r bin/
	$(RM) -r obj/

clean:cleanlocal




