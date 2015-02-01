CC0 = gcc
CC = $(CC0)
LD = $(CC0)
FLAGS1 = -O3
FLAGS2 = -lm
FLAGS3 = -std=c99 -Dlint -D__lint -Wall -Winline -Wno-deprecated -Wno-strict-overflow
DIR_H = ../share
FLAGS = $(FLAGS1) $(FLAGS3) -I$(DIR_H)

DEPS_SHARE = $(DIR_H)/ad_bucketio.h \
	$(DIR_H)/ad_fileio.h \
	$(DIR_H)/ad_lib.h \
	$(DIR_H)/ad_partition.h \
	$(DIR_H)/ad_print.h \
	$(DIR_H)/ad_random.h \
	$(DIR_H)/ad_readinput.h 

DEPS_SA1 = ad_lib_sa.h $(DIR_H)/ad_defs.h $(DEPS_SHARE)
DEPS_SA2 = $(DEPS_SA1)
DEPS_RSA = $(DEPS_SA1)

OBJS = $(patsubst %.h,%.o,$(DEPS_SHARE))
OBJS_SA1 = $(OBJS) ad_lib_sa.o
OBJS_SA2 = $(OBJS_SA1)
OBJS_RSA = $(OBJS_SA1)

all: ad_sa1 ad_sa2 ad_rsa

ad_lib_sa.o: ad_lib_sa.c ad_lib_sa.h $(DIR_H)/ad_lib.h $(DIR_H)/ad_bucketio.h $(DIR_H)/ad_defs.h
	$(CC) $(FLAGS) -c -o $@ $< $(FLAGS2)

ad_sa1: ad_sa1.c $(DEPS_SA1) $(OBJS_SA1)
	$(CC) $(FLAGS) -o $@.x $< $(OBJS_SA1) $(FLAGS2)

ad_sa2: ad_sa2.c $(DEPS_SA2) $(OBJS_SA2)
	$(CC) $(FLAGS) -o $@.x $< $(OBJS_SA2) $(FLAGS2)

ad_rsa: ad_rsa.c $(DEPS_RSA) $(OBJS_RSA)
	$(CC) $(FLAGS) -o $@.x $< $(OBJS_RSA) $(FLAGS2)

# Testing:
test:
	./utest.sh

# Cleaning:
clean c cl cle clea: 
	rm -f *.o *~ core *.x
cleano: 
	rm -f *.o *~ core

# End of file
