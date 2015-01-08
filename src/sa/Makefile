CC0 = gcc
CC = $(CC0)
LD = $(CC0)
FLAGS1 = -O3
FLAGS2 = -lm
FLAGS3 = -std=c99 -Dlint -D__lint -Wall -Winline -Wno-deprecated -Wno-strict-overflow
FLAGS = $(FLAGS1) $(FLAGS3)

DEPS_SHARE = ../share/ad_bucketio.h \
	../share/ad_fileio.h \
	../share/ad_lib.h \
	../share/ad_partition.h \
	../share/ad_print.h \
	../share/ad_random.h \
	../share/ad_readinput.h 

DEPS_SA1 = ad_lib_sa.h ../share/ad_defs.h $(DEPS_SHARE)
DEPS_SA2 = $(DEPS_SA1)
DEPS_RSA = $(DEPS_SA1)

OBJS = $(patsubst %.h,%.o,$(DEPS_SHARE))
OBJS_SA1 = $(OBJS) ad_lib_sa.o
OBJS_SA2 = $(OBJS_SA1)
OBJS_RSA = $(OBJS_SA1)

all: ad_sa1 ad_sa2 ad_rsa

ad_lib_sa.o: ad_lib_sa.c ad_lib_sa.h ../share/ad_lib.h ../share/ad_bucketio.h ../share/ad_defs.h
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
